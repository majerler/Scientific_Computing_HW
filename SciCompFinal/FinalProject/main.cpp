#include <iostream>
#include <cstdlib>
#include <ctime>
#include <autoencoder.h>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include<string>


#include <algorithm> // for std::copy

using namespace std;

vector<double> load_data(string file_name)
{
    vector<double> data;

    ifstream myfile;
    myfile.open ("scaled_data/" + file_name, ios::in);
    string line;
    while( getline(myfile, line) ) {
      double x = stod(line);
      data.push_back(x);
    }
    myfile.close();
    return data;
}


int main()
{
    AutoEncoder test(250, 50, 0.003);

    vector<vector<double>> vec;

    //Load Data
    for (int j = 0; j < 277; ++j)
    {
        vec.push_back( load_data("sample" + to_string(j) ) );
    }

    // K is the number of epochs
    for (int k = 0; k < 500; ++k) {
        double epoch_error = 0.0;
        // We go through every sample in our data set
        int samples = 277;
        for (int j= 0; j < samples; ++j)
        {
            double err = test.train( vec[j] );
            // We save the error from each sample
            epoch_error += err;
        }
        cout << "epoch "<< k << " average error " << epoch_error / (double) samples << endl;
    }

   // test.save();

    cout << "finshed!" << endl;

    return 0;

}
