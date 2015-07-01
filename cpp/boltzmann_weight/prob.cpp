#include "prob.h"

std::vector<double> Probability::Prob() {return prob;};

void Probability::ReadProb( std::string & filename )
{
	std::cout << "Reading " << filename << "..." << std::endl;
	std::string line;
	double value;
	double sum = 0;
	std::ifstream file (filename.c_str());
	if (file.is_open() )
	{
		while (file.good())
		{
			getline( file, line );
			if (not line.empty())
			{
				std::stringstream linestream( line );
				linestream >> value ;
				prob.push_back(value);
				sum += value;
			}
		}
	}
	else
	{
		std::cout << filename << " does not exist..." << std::endl;
		exit(1);
	}	
    // Normalizing the probability
	double inverseSum = 1.0/sum;
	for (int i=0; i<prob.size() ; i++)
	{
		prob[i] *= inverseSum;
	}
	return;
};

void Probability::clear()
{
    prob.clear();
}

std::vector<int> Bins::Bin() {return bin;}

void Bins::ReadBin( std::string & filename )
{
	std::cout << "Reading " << filename << "...\t";
	std::string line;
	int value;
    int nlines = 0;
	std::ifstream file (filename.c_str());
	if (file.is_open())
	{
		while (file.good())
		{
			getline( file, line);
			if (not line.empty()) 
			{
				std::stringstream linestream(line);
				linestream >> value;
				bin.push_back(value);
                nlines++;
			}
		}
	}
	else
	{
		std::cout << filename << " does not exist..." << std::endl;
		exit(1);
	}
    std::cout << "There are " << nlines << " points." << std::endl;
	return;
}

void Bins::clear()
{
    bin.clear();
};

std::vector< std::vector<double> > Dats::Dat(){return dat;}

void Dats::ReadDat( std::string & filename ) 
{
	std::cout << "Reading " << filename << "...\t"; 
	std::string line;
	std::vector<double> values;
	double n;
    int ncol;
    int nlines = 0;
	std::ifstream file (filename.c_str());
	if (file.is_open())
	{
		while (file.good())
		{
			getline( file, line);
			if ( not line.empty())
			{
                ncol = 0;
				std::stringstream linestream(line);
				while (linestream >> n) 
				{
                    std::stringstream value;
                    value >> n;
                    if (!value.good())
                    {
					    values.push_back(n);	
                        ncol++;
                    }
				}
                values.shrink_to_fit();
                if (values.size() > 0 )
                {
				    dat.push_back(values);
                    nlines++;
				    values.clear();
                }
                else 
                {
                    std::cout << "Not a number found on line " << nlines+1 << ".  Exiting..." << std::endl;
                    std::cout << line << std::endl;
                    exit(1);
                }
			}
		}
	}
	else
	{
		std::cout << filename << " does not exist..." << std::endl;
		exit(1);
	}
    std::cout << "There are " << ncol << " columns with " << nlines << " points." << std::endl;
	return;
};

void Dats::clear()
{
    dat.clear();
};

experiment Lists::ListDat(){return frame;}

void Lists::ReadList( std::string & filename )
{
    std::cout << "Reading " << filename << "." << std::endl;
    std::string line;
    std::vector<std::string> values; 
    std::string name;
    std::ifstream file(filename.c_str());
    if (file.is_open())
    {
        while (file.good())
        {
            getline( file, line);
            if (not line.empty())
            {
                std::stringstream linestream(line);
                while (linestream >> name)
                {
                    values.push_back(name);
                }
                list.push_back(values);
                values.clear();
            }
        }
    }
    else
    {
        std::cout << "Cannot open "<< filename << ". Exiting..." << std::endl;
        exit(1);
    }
    return;
};

void Lists::ReadListInputs()
{
    Bins bini;
    Dats dati;
    std::vector<int> bin;
    std::vector< std::vector<double> > dat;
    for (int i=0;i<list.size();i++)
    {
        
        bini.ReadBin( list[i][0] );
        dati.ReadDat( list[i][1] );
        bin = bini.Bin();
        dat = dati.Dat();
        if (dat.size() == bin.size()) 
        {
            for (int i=0; i< bin.size(); i++)
            {
                frame.bin.push_back(bin[i]);
                frame.dat.push_back(dat[i]);
            }
            bini.clear();
            dati.clear();
            bin.clear();
            dat.clear();
        }
        else
        {
            std::cout << "Sizes do not match! " << list[i][0] << "(" << bin.size() << ") " << list[i][1] << "(" << dat.size() << ")." << std::endl;
            // Should I make it still run?  That seems dangerously easy to get a bad output...
            exit(1);
        }
    }
    frame.dat.shrink_to_fit();
    frame.bin.shrink_to_fit();
    return;
};

void Lists::FrameProb(Probability probability)
{
    std::vector<double> prob = probability.Prob();
    int nbins = prob.size();
    std::cout << "\nThere are " << nbins << " bins." << std::endl;
    // Need to know how many times each bin is visited so that 
    // we can do frameprob = binprob/binvisits
    int * bincounts = new int[nbins];
    double * binprobs = new double[nbins];
    for (int i=0; i<nbins; i++)
    {
        bincounts[i] = 0;
        binprobs[i] = prob[i];
    }
    for (int i=0; i<frame.bin.size(); i++)
    {
        for (int j=0; j<nbins; j++)
        {
            if (frame.bin[i] == j)
            {
                bincounts[j] += 1;
                break;
            }
        }
    }
    // Here we do binprob/binvisits
    for (int i=0; i<nbins; i++)
    {
        if (bincounts[i] > 0)
        {
            binprobs[i] /= bincounts[i];
        }
    }
    delete [] bincounts;
    // Now get the probability of each frame
    for (int i=0; i<frame.bin.size(); i++)
    {
        for (int j=0; j<nbins; j++)
        {
            if (frame.bin[i] == j)
            {
                frame.prob.push_back(binprobs[j]);
                break;
            }
        }
    }
    delete [] binprobs;
    frame.prob.shrink_to_fit();
    return;
}

average Stats::angleAverage(experiment frame)
{
    double x, y, total, r, std, variance;
    double deg2rad = M_PI/180.;
    double rad2deg = 180./M_PI;
    // going to re-normalize the probability
    double probsum = 0.;
    for (int i=0; i< frame.prob.size(); i++)
    {
        probsum += frame.prob[i];
    }
    double inverseProbSum = 1./probsum;
    average result;
    for (int i=0; i<frame.dat[0].size(); i++)
    {
        x = 0;
        y = 0;
        for (int j=0; j<frame.dat.size(); j++)
        {
            x += cos(frame.dat[j][i]*deg2rad)*frame.prob[j]*inverseProbSum;
            y += sin(frame.dat[j][i]*deg2rad)*frame.prob[j]*inverseProbSum;
        }
        // Var = 1 - R
        // Var ranges from 0 to 1.
        // To make this analogous to angle, I should multiply by 360
        // rather than 180/pi, because this result is not a radian.
        // I knew the old variances looked too small...
        r = sqrt(x*x + y*y);
        variance = (1.-r)*360.;
        std = sqrt(variance);
        total = atan2(y,x)*rad2deg;
        result.avg.push_back(total);
        result.std.push_back(std);
    }
    return result;
};

average Stats::linearAverage(experiment frame)
{
    double total, vartotal, sum;
    // going to re-normalize the probability
    // This is mostly for bootstrapping
    double probsum = 0.;
    for (int i=0; i< frame.prob.size(); i++)
    {
        probsum += frame.prob[i];
    }
    double inverseProbSum = 1./probsum;
    average result;
    for (int i=0; i<frame.dat[0].size(); i++)
    {
        total = 0;
        vartotal = 0;
        sum = 0;
        for (int j=0; j<frame.dat.size(); j++)
        {
            total += frame.dat[j][i] * frame.prob[j] * inverseProbSum;
            sum += frame.prob[j];
        }
        for (int j=0; j<frame.dat.size(); j++)
        {
            vartotal += frame.prob[j] * inverseProbSum * frame.dat[j][i] * ( frame.dat[j][i] - total);
        }
        double std = sqrt(vartotal);
        result.avg.push_back(total);
        result.std.push_back(std);
    }
    return result;
};

void Stats::Average()
{
    boltzmannAverage.avg.clear();
    boltzmannAverage.std.clear();
    if (angular)
    {
        boltzmannAverage = angleAverage(frames);
    }
    else
    {
        boltzmannAverage = linearAverage(frames);
    }
    return;
}

average Stats::BoltzmannAverage(){return boltzmannAverage;}

average Stats::BootstrapAverage(){return bootstrapResults;}

void Stats::resampleAverage()
{
    average result;
    experiment resample;
    int randN;
    for (int i=0; i<frames.prob.size(); i++)
    {
        randN = rand() % frames.prob.size(); // generate random number between 0 and max number of frames
        resample.prob.push_back(frames.prob[randN]);
        resample.dat.push_back(frames.dat[randN]);
        resample.bin.push_back(frames.bin[randN]);
    }
    if (angular)
    {
        result = angleAverage(resample);
    }
    else
    {
        result = linearAverage(resample);
    }
    for (int i=0;i<result.avg.size();i++){
//        std::cout << "Resample: " << result.avg[i] << " +/- " << result.std[i] << std::endl;
    }
    results.push_back(result);

    resample.dat.clear();
    resample.prob.clear();
    resample.bin.clear();
    resample.dat.shrink_to_fit();
    resample.prob.shrink_to_fit();
    resample.bin.shrink_to_fit();
    return;
};

void Stats::bootstrap(int nresamples, bool useRandSeed)
{
    // Initiate random seed
    int seed = 47;
    if (useRandSeed)
    {
        seed = time(NULL);
    }
    std::cout << "Performing " << nresamples << " resampling experiments using seed: " << seed << std::endl;
    srand(seed);
    for (int i=0; i<nresamples; i++) 
    {
        resampleAverage();
    }
    // Let's make this an experiment so it can interact with our
    // previously made functions;
    experiment resampled;
    std::vector<double> datvec;
    double inverseNresamples = 1./nresamples;
    for (int i=0; i<results.size() ; i++){
        for (int j=0; j<results[0].avg.size(); j++){
            datvec.push_back(results[i].avg[j]);
        }
        resampled.dat.push_back(datvec);
        resampled.prob.push_back(inverseNresamples);
        datvec.clear();
    }

    if (angular)
    {
        bootstrapResults = angleAverage(resampled);
    }
    else
    {
        bootstrapResults = linearAverage(resampled);
    }
    for (int i=0;i<bootstrapResults.avg.size();i++){
    }
    std::cout<<std::endl;
}

void Stats::addFrame( experiment addframe, bool angle ){
    frames = addframe;
    angular = angle;
    Average();
}

void WriteOutputs::extension( std::string basename, std::string suffix, std::string & outname)
{
    std::string extension=".boltzmann";
    std::cout << basename.rfind(extension) << std::endl;
    if (basename.rfind(extension) != std::string::npos)
    {
        basename.replace(basename.rfind(extension),basename.length(),"");
    }
    outname.append(basename);
    outname.append(suffix);
    std::cout << outname << std::endl;
}

bool WriteOutputs::fexists( std::string filename )
{
	std::ifstream file (filename.c_str());
    if (file.good()) {
        file.close();
        return true;
    }
    else {
        file.close();
        return false;
    }
}

void WriteOutputs::backup( std::string filename )
{
    int n=0, nchars;
    std::string newname = filename;
    if (fexists(filename))
    {
        while (fexists(newname))
        {
            char buffer[1054];
            newname.clear();
            nchars = sprintf(buffer,"#%s.%d#",filename.c_str(),n);
            newname = buffer;
            n++;
            if (n>100)
            {
                std::cout << "Will not make more than 100 copies of " << filename << std::endl;
                std::cout << "Enter a different file name to backup to:" << std::endl;
                std::string filename;
                getline(std::cin,filename);
                newname = filename;
            }
        }
        std::cout << "Backing up " << filename << " to " << newname << std::endl;
        rename(filename.c_str(),newname.c_str());
    }
    std::cout << "Writing to " << filename << std::endl; 
    outfile = filename;
}

void WriteOutputs::write( std::string outfile, std::string probfile, average boltzmann, average bootstrap,int nsteps)
{
    backup(outfile);
    std::ofstream outputfile;
    outputfile.open(outfile.c_str());
    for (int i=0; i<boltzmann.avg.size() ; i++)
    {
        outputfile << i << " : " << probfile;
        outputfile << " avg/std ";
        outputfile << boltzmann.avg[i] << " ";
        outputfile << boltzmann.std[i];
        outputfile << "\n";
        
        std::cout << i << " : " << probfile;
        std::cout << " avg/std ";
        std::cout << boltzmann.avg[i] << " ";
        std::cout << boltzmann.std[i];
        std::cout << "\n";
    
        // This has problems with scientific notation
        /*
        outputfile << i << " : " << std::setw(30) << probfile;
        outputfile << std::setw(15) << " avg/std ";
        outputfile << std::setw(15) << std::setprecision(10) << boltzmann.avg[i];
        outputfile << std::setw(15) << std::setprecision(10) << boltzmann.std[i];
        outputfile << "\n";
        std::cout << i << " : " << std::setw(30) << probfile;
        std::cout << std::setw(15) << " avg/std ";
        std::cout << std::setw(15) << std::setprecision(10) << boltzmann.avg[i];
        std::cout << std::setw(15) << std::setprecision(10) << boltzmann.std[i];
        std::cout << "\n";
        */
    }
    int ndigits = log10(nsteps)+1;
    for (int i=0; i<bootstrap.avg.size() ; i++)
    {
        outputfile << i << " : " << "bootstrapping(" << nsteps << ") avg/ste ";
        outputfile << bootstrap.avg[i] << " ";
        outputfile << bootstrap.std[i]/sqrt(nsteps) << "\n";
        
        std::cout << i << " : " << "bootstrapping(" << nsteps << ") avg/ste ";
        std::cout << bootstrap.avg[i] << " ";
        std::cout << bootstrap.std[i]/sqrt(nsteps) << "\n";
        
        /*
        outputfile << i << " : " << std::setw(29-ndigits) << "bootstrapping(";
        outputfile << nsteps << ")";
        outputfile << std::setw(15) << " avg/ste ";
        outputfile << std::setw(15) << std::setprecision(10) << bootstrap.avg[i];
        outputfile << std::setw(15) << std::setprecision(10) << bootstrap.std[i]/sqrt(nsteps);
        outputfile << "\n";
        std::cout << i << " : " << std::setw(29-ndigits) << "bootstrapping(";
        std::cout << nsteps << ")";
        std::cout << std::setw(15) << " avg/ste ";
        std::cout << std::setw(15) << std::setprecision(10) << bootstrap.avg[i];
        std::cout << std::setw(15) << std::setprecision(10) << bootstrap.std[i]/sqrt(nsteps);
        std::cout << "\n";
        */
        /*
        std::cout << "\t99% likely for real value to be between ";
        std::cout << bootstrap.avg[i] - 2.5 * bootstrap.std[i]/sqrt(nsteps);
        std::cout << " and " << bootstrap.avg[i] + 2.5 * bootstrap.std[i]/sqrt(nsteps);
        std::cout << "\n";
         */
    }
    outputfile.close();
}

