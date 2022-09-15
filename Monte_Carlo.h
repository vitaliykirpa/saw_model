//
// Created by kamilla on 02.02.2022.
//

#ifndef SAW_MODELS_MONTE_CARLO_H
#define SAW_MODELS_MONTE_CARLO_H


#include "Model.h"

class Monte_Carlo {

public:
    Monte_Carlo (short int n, double J,
                 long int steps_to_equilibrium_, long int mc_steps_ ,
                 long int steps_to_write_);


    SAW_model *model;
    double J;
    short int L;
    long int steps_to_equilibrium;
    long int mc_steps;
    long int steps_to_write;

    virtual void run_simulation() = 0;

};

class Monte_Carlo_on_SAWS : public  Monte_Carlo
{
public:
    Monte_Carlo_on_SAWS (short int d, short int n, double J,
                 long int steps_to_equilibrium_ = 81000000, long int mc_steps_ = 5000000000000,
                 long int steps_to_write_=1000000);

    //~Monte_Carlo_on_SAWS();

    void run_simulation();

private:
    double p_for_local_update=0.5, p_for_reconnect=0.90;
};




#endif //SAW_MODELS_MONTE_CARLO_H
