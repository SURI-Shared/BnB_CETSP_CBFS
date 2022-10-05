#include <random>
#include <string>
#include <iostream>
#include <fstream>

using namespace std;

int main()
{
    string file_prefix = "CETSP-";
    vector<int> prob_sizes{ 50, 100, 150 };
    int prob_count = 10;
    
    
    for (int prob_size : prob_sizes)
    {
        for (int c = 1; c <= prob_count; ++c)
        {
			string filename = file_prefix + to_string(prob_size) + "-" + to_string(c) + ".txt";
            fstream ce;
			ce.open(filename, fstream::app);

			std::random_device rd;  // Will be used to obtain a seed for the random number engine
			std::mt19937 gen1(rd()); // Standard mersenne_twister_engine seeded with rd()
            std::mt19937 gen2(rd());
			std::uniform_real_distribution<> x(0.0, 16.0);
            std::uniform_real_distribution<> y(0.0, 10.0);

            for (int i = 0; i <= prob_size; ++i)
            {
                ce << x(gen1) << " " << y(gen2) << "\n";
            }

            ce.close();
        }
    }

    return 0;
}