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
// #include <FastNoise/Metadata.h>

 int benchmark(){
    auto node = FastNoise::New<FastNoise::FractalFBm>();
    std::cout << "SIMD level: " << node->GetSIMDLevel() << std::endl;
    auto noise = FastNoise::New<FastNoise::Perlin>();
    node->SetSource( noise );
    node->SetGain( 1.0 );
    node->SetOctaveCount(3);
    node->SetLacunarity(2.0);

    
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
        noise->GenUniformGrid2D( v.data(), 0, 0, N, N, 1./N, 1337 );

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
 
 enum Terrain { None=0, Sea, Swamp, Plain, Forest, Hills, Mountains};

 void save_BW_image(std::vector<float> &v, unsigned int N, std::string fn){
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
            int idx = N*y+x;
            float val = (v[idx]-min)*255/(max-min);
            image[4 * idx + 0] = val; // R
            image[4 * idx + 1] = val; // G
            image[4 * idx + 2] = val; // B
            image[4 * idx + 3] = 255; // alpha
        }
    }

    //write png
    unsigned error = lodepng::encode(fn, image, N, N);
    //if there's an error, display it
    if(error) std::cout << "lodePNG encoder error " << error << ": "<< lodepng_error_text(error) << "\t"<< fn << std::endl;
 }

 void save_Terrain_image(std::vector<Terrain> &v, unsigned int N, std::string fn){
    
    static unsigned char C_None[3] =   {0,   0,   0  };
    static unsigned char C_Sea[3] =   {31,  76,  191};
    static unsigned char C_Swamp[3] = {20,  115, 85 };
    static unsigned char C_Plain[3] = {185, 212, 68 };
    static unsigned char C_Forest[3]= {44,  110, 49 };
    static unsigned char C_Hills[3] = {148, 132, 92 };
    static unsigned char C_Mount[3] = {77 , 71,  61 };

    std::vector<unsigned char> image;
    image.resize(N * N * 4);
    for(unsigned y = 0; y < N; y++)
    {
        for(unsigned x = 0; x < N; x++) {
            int idx = N*y+x;
            unsigned char *color;
            switch(v[idx]){
                case Terrain::Sea:
                    color = C_Sea;
                    break;
                case Terrain::Swamp:
                    color = C_Swamp;
                    break;
                case Terrain::Plain:
                    color = C_Plain;
                    break;
                case Terrain::Forest:
                    color = C_Forest;
                    break;
                case Terrain::Hills:
                    color = C_Hills;
                    break;
                case Terrain::Mountains:
                    color = C_Mount;
                    break;
                default:
                    color = C_None;
                    break;
            }
            image[4 * idx + 0] = color[0]; // R
            image[4 * idx + 1] = color[1]; // G
            image[4 * idx + 2] = color[2]; // B
            image[4 * idx + 3] = 255; // alpha
        }
    }

    //write png
    unsigned error = lodepng::encode(fn, image, N, N);
    //if there's an error, display it
    if(error) std::cout << "lodePNG encoder error " << error << ": "<< lodepng_error_text(error) << "\t"<< fn << std::endl;
 }



 int main(){

    int seed_sea = 1253;            // seas
    int seed_mnt = seed_sea+51;     // mountains
    int seed_per = seed_sea+145;    // percipitation

    float freq_sea = 0.003;
    float freq_mnt = freq_sea * 0.7;
    float freq_per = freq_sea * 2.5;

    auto noise = FastNoise::New<FastNoise::Perlin>();

    // seas
    auto seas = FastNoise::New<FastNoise::FractalFBm>();
    seas->SetSource( noise );
    seas->SetGain( 1.0 );
    seas->SetOctaveCount(3);
    seas->SetLacunarity(2.0);

    // mountains
    auto  mnts = FastNoise::New<FastNoise::FractalRidged>();
    mnts->SetSource( noise );
    mnts->SetGain( 1.0 );
    mnts->SetOctaveCount(3);
    mnts->SetLacunarity(2.2);

    // percipitation
    auto  pers = FastNoise::New<FastNoise::FractalFBm>();
    pers->SetSource( noise );
    pers->SetGain( 1.0 );
    pers->SetOctaveCount(3);
    pers->SetLacunarity(3.0);

    
    unsigned int N = 512;
    std::vector<float> map_seas, map_mnts, map_pers;
    map_seas.resize(N * N);
    map_mnts.resize(N * N);
    map_pers.resize(N * N);

    std::vector<Terrain> map_out;
    map_out.resize(N * N);

    // noise generation
    std::clock_t c_start = std::clock();
    seas->GenUniformGrid2D( map_seas.data(), 0, 0, N, N, freq_sea, seed_sea );
    mnts->GenUniformGrid2D( map_mnts.data(), 0, 0, N, N, freq_mnt, seed_mnt );
    pers->GenUniformGrid2D( map_pers.data(), 0, 0, N, N, freq_per, seed_per );
    std::clock_t c_end = std::clock();
    float time_elapsed_ms = 1000.0 * (c_end-c_start) / CLOCKS_PER_SEC;
    std::cout << "Noise generation:  CPU time used: " << time_elapsed_ms << " ms for "<< N*N <<" points\n";


    
    c_start = std::clock();

    // adjust seas to make an island
    float center = N*0.5;
    float sea_edge_min = N*0.3 * N*0.3; // keep these squared for comparison with rsq
    float sea_edge_max = N*0.45 * N*0.45;
    float sea_edge_width = sea_edge_max - sea_edge_min;
    for(unsigned y = 0; y < N; y++)
    {
        float ydif = y-center;
        ydif*=ydif;
        for(unsigned x = 0; x < N; x++) {
            float xdif = x-center;
            xdif*=xdif;
            float rsq = ydif + xdif;
            if(rsq>sea_edge_min){
                float shift =(rsq - sea_edge_min)/sea_edge_width;
                // float scale = 1.0 - 0.7*shift;
                // map_seas[N*y+x]*=scale;
                map_seas[N*y+x]-= 0.5*shift;
            }
        }
    }

    // terrain assignment
    for(unsigned i = 0; i < N*N; ++i){
        
        if(map_seas[i]<-0.2){
            map_out[i]=Terrain::Sea;
        }
        else if(map_mnts[i]>0.9 && map_seas[i]>0){
            map_out[i]=Terrain::Mountains;
        }
        else if(map_mnts[i]>0.75 && map_seas[i]>0){
            map_out[i]=Terrain::Hills;
        }
        else if(map_pers[i]>0.55 && map_mnts[i] < 0.2 ){
            map_out[i]=Terrain::Swamp;
        }
        else if(map_pers[i]>0.15){
            map_out[i]=Terrain::Forest;
        }
        else{
            map_out[i]=Terrain::Plain;
        }
    }
    c_end = std::clock();
    time_elapsed_ms = 1000.0 * (c_end-c_start) / CLOCKS_PER_SEC;
    std::cout << "Terrain Assignment:  CPU time used: " << time_elapsed_ms << " ms for "<< N*N <<" points\n";


    // save images
    std::ostringstream seas_fn;
    seas_fn << "test_seas_" << N << ".png";
    save_BW_image(map_seas, N, seas_fn.str());

    std::ostringstream mnts_fn;
    mnts_fn << "test_mnts_" << N << ".png";
    save_BW_image(map_mnts, N, mnts_fn.str());

    std::ostringstream pers_fn;
    pers_fn << "test_pers_" << N << ".png";
    save_BW_image(map_pers, N, pers_fn.str());

    std::ostringstream out_fn;
    out_fn << "test_out_" << N << ".png";
    save_Terrain_image(map_out, N, out_fn.str());

    
    
    
    return(0);
 }
 