#ifndef AUTOENCODER_H
#define AUTOENCODER_H

#include <vector>
#include <iostream>
#include <net_tools.h>
#include <fstream>
#include <string.h>

#include <random>
//#include <cstdlib>
//#include <ctime>


class AutoEncoder
{
private:
    int inputDim;
    int hiddenDim;
    double learningRate;

    std::vector<double> W_hidden;
    std::vector<double> W_out;
    std::vector<double> b_hidden;
    std::vector<double> b_out;
    std::vector<double> error_hidden;
    std::vector<double> error_out;
    std::vector<double> z_hidden;
    std::vector<double> z_out;
    std::vector<double> a;
    std::vector<double> y_hat;

    void forward(  std::vector<double> sample);
    void backward( std::vector<double> samp);
    void update(   std::vector<double> samp);

public:
    AutoEncoder();
    AutoEncoder(int inputDim_, int hiddenDim_, double learningRate_);
    double train(std::vector<double> samp);
    void print();
    void save();
};

#endif // AUTOENCODER_H
