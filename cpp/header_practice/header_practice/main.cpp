//
//  main.cpp
//  header_practice
//
//  Created by Andrew Ritchie on 4/5/13.
//  Copyright (c) 2013 Andrew Ritchie. All rights reserved.
//

#include <iostream>
#include "File.h"

using namespace std;

int main(int argc, const char * argv[])
{

    myfunction(10);
    MyClass awesome(8);
    std::cout << awesome.your_number() << std::endl;
    awesome.class_function(15);
    std::cout << awesome.your_number() << std::endl;

    return 0;
}

