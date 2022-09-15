//
// Created by kamilla on 01.02.2022.
//

#ifndef SAW_MODELS_MODEL_H
#define SAW_MODELS_MODEL_H

#include <random>
#include <map>
#include "observable.h"
#include "Lattice.h"

const double PI = std::atan(1.0)*4;


class Lattice;

class Model{
public:
    virtual long int number_of_spins() = 0;
    double Energy_value() {return E;};
    std::valarray<double> sequence_on_lattice;
    /*template<typename T>
     T Hamiltonian ();*/ //no idea for double/int energy in inherited
    virtual void save_measurements() = 0 ;
    virtual void write_file(long int i, double J) = 0;

    //virtual void save_counts() = 0;
protected:
    long int L;
    Lattice *lattice;
    double E;
    virtual void Energy () = 0; //call Hamiltonian here

    mc_stats::ScalarObservable<double> energy; //сохранение энергии
    mc_stats::ScalarObservable<double> energy_sq;
    mc_stats::ScalarObservable<double> energy_4;

    double h_l = 0.001; //длина  бина

protected:
    std::map <long int, long int> count_E;
    std::map <long int, long int> count_m2;

};

class SAW_model : public Model
{
public:
    SAW_model(short int n, short d);
    int check_dimenstion() {return dimension;};
    std::valarray<int> directions; //их n-1;
    std::valarray<long int> next_monomers;
    std::valarray<long int> previous_monomers;
    long int end_conformation=0, start_conformation=0;
    virtual void Energy() = 0;
    virtual long int number_of_spins() = 0;

    void Reconnect(int j);
    virtual void FlipMove (double J) = 0;
    virtual void ClusterStep (double J) = 0;

    void save_geometry();
    virtual void save_measurements() = 0;
    void calc_bulc();
    //virtual void save_counts() ;

    mc_stats::ScalarObservable<double> dists;
    mc_stats::ScalarObservable<double> gyration;
    //можно ли для 3D наблюдать за плоскостью?
    mc_stats::ScalarObservable<double> eigval1;
    mc_stats::ScalarObservable<double> eigval2;
    mc_stats::ScalarObservable<double> aratio;

    //добавить заранее балки для 3Д?
    mc_stats::ScalarObservable<double> bulk6;
    mc_stats::ScalarObservable<double> bulk5;
    mc_stats::ScalarObservable<double> bulk4; //доля узлов с 4 соседами
    mc_stats::ScalarObservable<double> bulk3;
    mc_stats::ScalarObservable<double> bulk2;

    std::map <long int, long int> count_R2;
    std::map <long int, long int> count_X;
    std::map <long int, long int> count_Y;

protected:
    int dimension;


};


class XY_SAW : public SAW_model{
public:
    XY_SAW (short int n, short d = 2) ;
    long int number_of_spins() { return L;};
    void Energy ();

    void FlipMove (double J);
    void ClusterStep (double J);

    void save_measurements();
    //void save_counts();
    void write_file(long int i, double J);

protected:
    mc_stats::ScalarObservable<double> mags_sin;
    mc_stats::ScalarObservable<double> mags_cos;
    mc_stats::ScalarObservable<double> magnetization_sq;
    mc_stats::ScalarObservable<double> magnetization_4;
    std::map <long int, long int> count_cos;

    std::valarray<bool> used_coords;

};


#endif //SAW_MODELS_MODEL_H
