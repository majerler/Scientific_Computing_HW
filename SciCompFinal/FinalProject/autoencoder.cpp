#include "autoencoder.h"


AutoEncoder::AutoEncoder()
{

}

AutoEncoder::AutoEncoder(int inputDim_, int hiddenDim_, double learningRate_)
{
    inputDim = inputDim_ ;
    hiddenDim = hiddenDim_ ;
    learningRate = learningRate_ ;

    // We want W_hidden and W_out to be matrices of size (input*hidden) and (hidden*input)
    W_hidden.resize(inputDim*hiddenDim);
    W_out.resize(hiddenDim*inputDim);

    // We will add the bias to the multiplication so we want dims to match
    b_hidden.resize(hiddenDim);
    b_out.resize(inputDim);

    // We need a place to sstore the errors we will calculate
    error_hidden.resize(hiddenDim);
    error_out.resize(inputDim);

    z_hidden.resize(hiddenDim);
    z_out.resize(inputDim);

    a.resize(hiddenDim);
    y_hat.resize(inputDim);


    std::default_random_engine generator (42);
    std::normal_distribution<double> distribution (0.0, 1.0);

//    for (int j = 0; j < 5; ++j) {
//        std::cout << (double) distribution(generator) << std::endl;
//    }

    for (int j = 0; j < inputDim*hiddenDim; ++j) {
        W_hidden[j] = (double) distribution(generator)  ;
        W_out[j] =    (double) distribution(generator)  ;
    }

    for (int j = 0; j < hiddenDim; ++j) {
        b_hidden[j] = (double) distribution(generator)  ;
    }

    for (int j = 0; j < inputDim; ++j) {
        b_out[j] =    (double) distribution(generator)  ;
    }

}

double AutoEncoder::train(std::vector<double> samp)
{
    forward(samp);
    backward(samp);
    update(samp);

    // Calculate the MSE
    double error = 0.0;
    for (int j = 0; j < y_hat.size(); ++j) {
        error += 0.5 * ( y_hat[j] - samp[j] ) * ( y_hat[j] - samp[j] );
    }

    return error / (double) y_hat.size();

}

void AutoEncoder::forward(std::vector<double> sample)
{

    // Reset values
        for (int j = 0; j < hiddenDim ; ++j) {
            z_hidden[j] = 0.  ;
            a[j]        = 0. ;
        }

        for (int j = 0; j < inputDim ; ++j) {
            z_out[j] = 0. ;
            y_hat[j] = 0. ;
        }

    // Multiply input by hidden matrix
    for (int i= 0; i < hiddenDim; ++i) {
        for (int j = 0; j < inputDim ; ++j) {
            z_hidden[i] += W_hidden[j + i*inputDim] * sample[j] ;
        }
    }

    // Add bias and use activation function
    for (int j = 0; j < hiddenDim ; ++j) {
        z_hidden[j] = z_hidden[j] + b_hidden[j];
        a[j] = relu( z_hidden[j] );
    }

    // multiply by output matrix
    for (int i= 0; i < inputDim; ++i) {
        for (int j = 0; j < hiddenDim ; ++j) {
            z_out[i] +=  W_out[j + i*hiddenDim] * a[j];
        }
    }

    // Add bias and use activation function
    for (int i = 0; i < inputDim ; ++i) {
        z_out[i] = z_out[i] + b_out[i];
        y_hat[i] = sigmoid( z_out[i] );
    }

}

void AutoEncoder::backward(std::vector<double> samp)
{
    // reset values for error_out and error_hidden!
    for (int i = 0; i < inputDim; ++i) {
        error_out[i] = 0.  ;
    }

    for (int i = 0; i < hiddenDim; ++i) {
        error_hidden[i] = 0. ;
    }

    // calculate error out
    for (int i = 0; i < inputDim; ++i) {
        error_out[i] = (y_hat[i] - samp[i]) * sigmoid_dx( z_out[i] );
    }

    // Calculate hidden error!
    for (int i= 0; i < hiddenDim; ++i) {
        for (int j = 0; j < inputDim ; ++j) {
            error_hidden[i] += error_out[j]  * W_out[i + j*hiddenDim] ;
        }
    }
    for (int j = 0; j < hiddenDim ; ++j) {
        error_hidden[j] = error_hidden[j] * relu_dx( z_hidden[j] );
    }

}

void AutoEncoder::update(std::vector<double> samp)
{
    //Update hidden weights
    for (int i= 0; i < hiddenDim; ++i) {
        for (int j = 0; j < inputDim ; ++j) {
            W_hidden[j + i*inputDim] -= samp[j] * error_hidden[i] * learningRate ;
        }
    }

    //Update hidden bias
    for (int i= 0; i < hiddenDim; ++i) {
        b_hidden[i] -= error_hidden[i] * learningRate ;
    }

    // Update Output Weights
    for (int i= 0; i < inputDim; ++i) {
        for (int j = 0; j < hiddenDim ; ++j) {
            W_out[j + i*hiddenDim] -= a[j] * error_out[i] * learningRate;
        }
    }

    // Update Output Weights
    for (int i= 0; i < inputDim; ++i) {
        b_out[i] -=  error_out[i] * learningRate;
    }
}


void AutoEncoder::print()
{
    int N = inputDim ;
    int M = hiddenDim;

    std::cout << "Weight Matrix - Hidden" << std::endl;
    for (int q = 0; q < M; ++q) {
        for (int p = 0; p < N; ++p) {
            std::cout << W_hidden[ p + q * N] << " ";
        }
         std::cout << std::endl;
    }

    std::cout << std::endl;

    std::cout << "Bias - Hidden" << std::endl;
    for (int q = 0; q < M; ++q) {
            std::cout << b_hidden[q] << " ";
    }

    std::cout << std::endl;

    std::cout << "Weight Matrix - Output" << std::endl;
    for (int q = 0; q < N; ++q) {
        for (int p = 0; p < M; ++p) {
            std::cout << W_out[ p + q * M] << " ";
        }
         std::cout << std::endl;
    }

    std::cout << std::endl;

    std::cout << "Bias - Output" << std::endl;
    for (int q = 0; q < N; ++q) {
            std::cout << b_out[ q ] << " ";
    }

}

void AutoEncoder::save()
{
    int N = inputDim ;
    int M = hiddenDim;

    std::ofstream w_hid;
    w_hid.open ("W_hidden.txt", std::ios::out | std::ios::app);

    for (int a = 0; a < N*M; ++a) {
        w_hid << std::to_string(W_hidden[a]) + "\n";
    }

    w_hid.close();

    std::ofstream w_out;
    w_out.open ("W_output.txt", std::ios::out | std::ios::app);

    for (int a = 0; a < N*M; ++a) {
        w_out << std::to_string(W_out[a]) + "\n";
    }

    w_out.close();

    std::ofstream b_hid;
    b_hid.open ("b_hidden.txt", std::ios::out | std::ios::app);

    for (int a = 0; a < M; ++a) {
        b_hid << std::to_string(b_hidden[a]) + "\n";
    }

    b_hid.close();

    std::ofstream b_o;
    b_o.open ("b_output.txt", std::ios::out | std::ios::app);

    for (int a = 0; a < N; ++a) {
        b_o << std::to_string(b_out[a]) + "\n";
    }

    b_o.close();

}

