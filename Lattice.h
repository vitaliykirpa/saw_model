//
// Created by kamilla on 01.02.2022.
//

#include <valarray>
#include <vector>
#include "Model.h"

#ifndef SAW_MODELS_LATTICE_H
#define SAW_MODELS_LATTICE_H

typedef unsigned  long int coord_t;

class Model;
class SAW_model;

class Lattice {
public:
    Lattice( short max_seq_size = 10 ) {lattice_side = max_seq_size; };
    long int lattice_size() {return lattice_side;};
    long int NumberOfNodes () {return number_of_nodes;};
    inline const virtual short int ndim() = 0 ;
    inline const virtual short int ndim2() = 0;
    std::valarray<coord_t> map_of_contacts_int;
    //that is a mistake
    //solve it later
    std::valarray<int> inverse_steps = {1,0,3,2} ;//{ 1, 0, 3, 2, 5, 4 }; // = {1,0,3,2} ;//= { 1, 0, 3, 2, 5, 4 };//= {1,0,3,2} ;
    //std::vector<int>  inverse_steps ;//= std::vector<int>(10);//? should work for different lattices
    //std::valarray<int> inverse_steps ;
    std::vector <std::vector<int>> steps ;
    virtual void radius(SAW_model  *model) = 0;
protected:
    short lattice_side = 0;
    int number_of_nodes = 0;
    virtual void create_lattice() = 0;
};

class Lattice_2D : public Lattice {
public:
    Lattice_2D( short max_seq_size = 10);
    inline const short int ndim() { return 2; }
    inline const short int ndim2() {return 4;}

    //std::valarray<int> inverse_steps = {1,0,3,2} ;
    std::vector <std::vector<int>> steps = {{1,0}, {-1,0}, {0,1}, {0,-1} };
    void radius(SAW_model  *model);
private:
    void create_lattice();
};


class Lattice_3D : public Lattice {
public:
    Lattice_3D( short max_seq_size = 10);
    inline const short int ndim() { return 3; }
    inline const short int ndim2() {return 6;}

    //std::valarray<int> inverse_steps = {1,0,3,2} ;
    std::vector <std::vector<int>> steps {   {1,0,0}, {0,1,0}, {0,0,1},
                                               {0,0,-1}, {0,-1,0}, {-1,0,0} };
    void radius(SAW_model  *model);
private:
    void create_lattice();
};

#endif //SAW_MODELS_LATTICE_H
