#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <ctime>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>


#include "lodepng.h"
#include "FastNoise/FastNoise.h"


 int main(){


    auto node = FastNoise::New<FastNoise::FractalFBm>();
    std::cout << node->GetSIMDLevel() << std::endl;
    auto noise = FastNoise::New<FastNoise::Perlin>();
    node->SetSource( noise );

    
    unsigned int nums[6]={16, 32, 64, 128, 256, 512};
    // float point[2];
    // unsigned int nums[1]={512};
    std::ofstream benchlog;
    benchlog.open ("bench.log");
    for(unsigned int benchloop=0; benchloop<6; ++benchloop){
    // for(unsigned int benchloop=0; benchloop<1; ++benchloop){
        unsigned int N = nums[benchloop];
        std::vector<float> v;
        v.resize(N * N);
        std::clock_t c_start = std::clock();
        // for(unsigned y = 0; y < N; y++)
        // {
        //     point[1] = y*8./N;
        //     for(unsigned x = 0; x < N; x++) {
        //         point[0] = x*8./N;
        //         // v[N*y+x] = test_layer.getValue(point);
        //     }
        // }

        node->GenUniformGrid2D( v.data(), 0, 0, N, N, N/8., 1337 );

        std::clock_t c_end = std::clock();
        float time_elapsed_ms = 1000.0 * (c_end-c_start) / CLOCKS_PER_SEC;
        std::cout << "2D periodic:  CPU time used: " << time_elapsed_ms << " ms for "<< N*N <<" points\n";
        benchlog << N*N << "\t" << time_elapsed_ms << std::endl;

        float min=1;
        float max=-1;
        for(unsigned i = 0; i < N*N; ++i)
        {
            if(v[i]>max) max=v[i];
            if(v[i]<min) min=v[i];
        }
        std::cout <<"min: " << min << "\tmax: "<< max << std::endl;

        std::vector<unsigned char> image;
        image.resize(N * N * 4);
        for(unsigned y = 0; y < N; y++)
        {
            for(unsigned x = 0; x < N; x++) {
                float val = (v[N*y+x]-min)*255/(max-min);
                image[4 * N * y + 4 * x + 0] = val; // R
                image[4 * N * y + 4 * x + 1] = val; // G
                image[4 * N * y + 4 * x + 2] = val; // B
                image[4 * N * y + 4 * x + 3] = 255; // alpha
            }
        }

        //write png
        std::ostringstream out;
        out << "test_2D_" << N << "x" << N << ".png";
        unsigned error = lodepng::encode(out.str(), image, N, N);
        //if there's an error, display it
        if(error) std::cout << "lodePNG encoder error " << error << ": "<< lodepng_error_text(error) << std::endl;

    
    }
    benchlog.close();

    return(0);
 }
 