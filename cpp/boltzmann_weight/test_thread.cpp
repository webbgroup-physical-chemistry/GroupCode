// thread example
#include <iostream>       // std::cout
#include <thread>         // std::thread
#include <vector>

struct F
{
    void operator() () const
    {
        std::cout<<"Printing from another thread "<<std::endl;
    }
};

void func(int i,int it)
{
    std::cout << "\n" << i << ") on THREAD " << it ;
}
void func2(int i,int it)
{
 //   std::cout << "\n" << i << ") on THREAD " << it ;
    printf("\n%d on THREAD %d\n",i,it);
    for (int j=0;j<1e8;j++){i+=j;}
}

int main()
{
    int nthreads = std::thread::hardware_concurrency();
    std::vector<std::thread> threads;
/*    for (int i=0;i<nthreads;i++){
        threads.push_back(std::thread(func2,0,i));
    }
    for (auto& th : threads) th.join();
*/
    int n=0;
    for (int i=0;i<nthreads*10;i+=nthreads){
        for (int j=0;j<nthreads;j++){
//            std::cout << i << " " << n << " " << j << std::endl;
            threads.push_back(std::thread(func2,n,j));
            n++;
        }
        for (auto& th : threads) th.join();
        threads.clear();

    }
    /*
    for (int j=0; j<nthreads; j++){
        for (int i=0;i<nthreads*5;i++){
            threads.at(j) = std::thread(func2,i,j);
        }
        threads.at(j).join();
        */
//        std::thread t(func,i,1);
//        std::thread h(func2,i,2);
//        t.join();
//        h.join();

    return 0;
}
