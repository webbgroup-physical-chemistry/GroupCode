//
//  File.cpp
//  header_practice
//
//  Created by Andrew Ritchie on 4/5/13.
//  Copyright (c) 2013 Andrew Ritchie. All rights reserved.
//

#include "File.h"

void myfunction(int number)
{
    std::cout << "Hello World! You entered " << number << std::endl;
}

MyClass::MyClass(int number)
{
    my_number=number;
}
int MyClass::your_number() {return my_number;}
void MyClass::class_function(int number)
    {
        std::cout << "This is a class.  You entered the number " << number << std::endl;
        my_number=number;
    }
