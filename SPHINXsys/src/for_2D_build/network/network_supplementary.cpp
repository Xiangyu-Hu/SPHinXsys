/**
 * @file 	branch.cpp
 * @brief 	Here, Functions belong to BaseBranch are given.
 * @author	chi ZHang and Xiangyu Hu
 * @version	0.2.0
 */
 #include "network.h"
 //=================================================================================================//
namespace SPH
{
    //=================================================================================================//
    Vecd Node::getGradient(Point pt, Real delta)
    {

        std::cout << "Only work for 3 dimension !" << std::endl;
        std::cout << __FILE__ << ':' << __LINE__ << std::endl;
		exit(1);
        return Vecd(0);
    }
    //=================================================================================================//
}
//=================================================================================================//