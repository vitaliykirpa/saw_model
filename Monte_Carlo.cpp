//
// Created by kamilla on 02.02.2022.
//

#include "Monte_Carlo.h"

Monte_Carlo::Monte_Carlo(short n_, double J_, long steps_to_equilibrium_, long mc_steps_, long steps_to_write_)
{
    L = n_;
    J = J_;
    steps_to_equilibrium = steps_to_equilibrium_;
    mc_steps = mc_steps_;
    steps_to_write = steps_to_write_;
}

Monte_Carlo_on_SAWS::Monte_Carlo_on_SAWS(short d, short n_, double J_, long steps_to_equilibrium_, long mc_steps_, long steps_to_write_)
                                         : Monte_Carlo(n_,J_,steps_to_equilibrium_,mc_steps_,steps_to_equilibrium_)
{
    model = new XY_SAW (n_, d);

}

void Monte_Carlo_on_SAWS::run_simulation()
{

    std::random_device generator;
    std::uniform_real_distribution<double> distribution(0.0,1.0);
//std::mt19937 generator(1234);
    std::random_device generators1;
    std::uniform_int_distribution<int> distribution1(0, 3);
//std::mt19937 generators1(123);


    double typeOfUpdate; //0 - простой; 1 - реконнект
    int step_on_lattice;//выбор одного из соседей


    long int all_steps=steps_to_equilibrium+mc_steps;

    for (long int i=0; i<all_steps+2; i++) {
        //std::cout << E << std::endl;
        //std::cout << "STEP : " << i << std::endl;
        typeOfUpdate = distribution(generator);
        if (typeOfUpdate < p_for_local_update) {
            model->FlipMove(J);

        }
        else if (typeOfUpdate<p_for_reconnect) {

            step_on_lattice = distribution1(generators1);
            model->Reconnect(step_on_lattice);

        }
        else
        {
            model->ClusterStep(J);
        }
        if (  i > steps_to_equilibrium &&  i%steps_to_write==0    )
        {
            model->save_geometry();
            model->calc_bulc();
            model->save_measurements();
        }

        if ( i> steps_to_equilibrium && i%(100*steps_to_write)==0 )
        {
            model->write_file(i, J);
        }
    }
}

/*Monte_Carlo_on_SAWS::~Monte_Carlo_on_SAWS()
{
    delete model;

}*/


//#include "gperftools/profiler.h"
//unsigned int thread_qty = std::max(atoi(std::getenv("OMP_NUM_THREADS")), 1);
//omp_set_num_threads(1);
//#include <omp.h>

int main(int argc, char *argv[])
{

    //omp_set_num_threads(1);
    //Lattice_2D l(100);
    int N = std::stoi(argv[1]);
    double J = 0.001*(double)std::stoi(argv[2]);


    //Monte_Carlo_on_SAWS m(3,N,J,10000,100000000,100);
    Monte_Carlo_on_SAWS m(2,N,J,10000,100000000000000,100000);
    //ProfilerStart("main.prof");
    m.run_simulation();
    //ProfilerStop();
    //XY_SAW xy1(10);

    return 0;
}


