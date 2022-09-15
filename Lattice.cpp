//
// Created by kamilla on 01.02.2022.
//

#include "Lattice.h"



Lattice_2D::Lattice_2D(short max_seq_size) : Lattice(max_seq_size) {
    //lattice_side = max_seq_size;
    number_of_nodes = lattice_side*lattice_side;
    create_lattice();
}

Lattice_3D::Lattice_3D(short max_seq_size) : Lattice(max_seq_size) {
    //lattice_side = max_seq_size;
    number_of_nodes = lattice_side*lattice_side*lattice_side;
    create_lattice();
}

void Lattice_2D::create_lattice() {

    long int x, y;
    div_t n;
    map_of_contacts_int.resize(lattice_side*lattice_side*ndim2());
    for (long int i =0; i<number_of_nodes ; i++){
        map_of_contacts_int[ndim2()*i] = i+1;
        map_of_contacts_int[ndim2()*i+1] = i-1;
        map_of_contacts_int[ndim2()*i+2] = i+lattice_side;
        map_of_contacts_int[ndim2()*i+3] = i-lattice_side;
        n=div(i, lattice_side);
        x=n.rem;
        y=n.quot;
        for (int j =0; j<ndim2(); j++){
            if(x==0){
                map_of_contacts_int[ndim2()*i+1] = i+lattice_side-1;
            }
            if(x==(lattice_side-1)){
                map_of_contacts_int[ndim2()*i] = i-(lattice_side-1);
            }
            if(y==0){
                map_of_contacts_int[ndim2()*i+3] = lattice_side*(lattice_side-1)+x;
            }
            if(y==(lattice_side-1)){
                map_of_contacts_int[ndim2()*i+2] = x;
            }
        }
    }

}

void Lattice_3D::create_lattice() {
    coord_t x, y,z;
    div_t n;
    map_of_contacts_int.resize(lattice_side*lattice_side*lattice_side*ndim2());
    coord_t l;
    for (long int i =0; i<number_of_nodes ; i++){
        map_of_contacts_int[ndim2() * i] = i + 1;
        map_of_contacts_int[ndim2() * i + 1] = i - 1;
        map_of_contacts_int[ndim2() * i + 2] = i + lattice_side;
        map_of_contacts_int[ndim2() * i + 3] = i - lattice_side;

        map_of_contacts_int[ndim2() * i + 4] = i + lattice_side * lattice_side;
        map_of_contacts_int[ndim2() * i + 5] = i - lattice_side * lattice_side;

        l = lattice_side * lattice_side;
        n = div(i, l);
        z = n.quot;
        n = div( n.rem, lattice_side);
        x = n.rem;
        y = n.quot;
        //x = n.rem % lattice_side; //n.rem;
        //y = n.rem / lattice_side; //n.quot;
        for (int j = 0; j < ndim2(); j++) {
            if (x == 0) {
                map_of_contacts_int[6 * i + 1] = i + lattice_side - 1;
            }
            if (x == (lattice_side - 1)) {
                map_of_contacts_int[6 * i] = i - (lattice_side - 1);
            }
            if (y == 0) {
                //map_of_contacts_int[6 * i + 3] = lattice_side * (lattice_side - 1) + x;
                map_of_contacts_int[6 * i + 3] =  lattice_side * (lattice_side - 1) + i;
            }
            if (y == (lattice_side - 1)) {
                //map_of_contacts_int[6 * i + 2] = x;
                map_of_contacts_int[6 * i + 2] = i - lattice_side * (lattice_side - 1) ;
            }

            if (z == 0) {
                //map_of_contacts_int[6 * i + 5] = x + lattice_side * y + lattice_side * lattice_side * (lattice_side - 1);
                map_of_contacts_int[6 * i + 5] = i + lattice_side * lattice_side * (lattice_side - 1);
            }
            if (z == lattice_side - 1) {
                //map_of_contacts_int[6 * i + 4] = x + lattice_side * y;
                map_of_contacts_int[6 * i + 4] = i - lattice_side * lattice_side * (lattice_side - 1);
            }
        }
    }
}

void Lattice_3D::radius(SAW_model *model) {
    long int point1x = model->end_conformation % lattice_side;
    long int point1z = model->end_conformation / (lattice_side * lattice_side);
    long int point1y = (model->end_conformation % (lattice_side * lattice_side)) / lattice_side;
    long int point1xs = model->start_conformation % lattice_side;
    long int point1zs = model->start_conformation / (lattice_side * lattice_side);
    long int point1ys = (model->start_conformation % (lattice_side * lattice_side)) / lattice_side;

    //расстояние на торе
    long int xdiff = abs(point1x- point1xs);
    if (xdiff > (lattice_side  / 2))
        xdiff = lattice_side - xdiff;

    long int ydiff = abs(point1y- point1ys);
    if (ydiff > (lattice_side / 2))
        ydiff = lattice_side - ydiff;

    long int zdiff = abs(point1z- point1zs);
    if (zdiff > (lattice_side / 2))
        zdiff = lattice_side - zdiff;

    long int r = xdiff *xdiff  + ydiff*ydiff + zdiff*zdiff;
    model->dists << r;

    model->count_R2[r]+=1;
    //model->count_X[xdiff]+=1;
    //model->count_Y[ydiff]+=1;

    std::vector <long long int> xs,ys,zs;
    xs.push_back(0);
    ys.push_back(0);
    zs.push_back(0);
    int direction;


    long int current = model->start_conformation;

    //std::vector <std::vector<int>> steps = {   {1,0}, {-1,0}, {0,1}, {0,-1} };
    for (long long int i =1; i < model->number_of_spins(); i++)
    {
        direction = model->directions[current];
        long long int x = xs.back()+steps[direction][0];
        long long int y = ys.back()+steps[direction][1];
        long long int z = zs.back()+steps[direction][2];
        xs.push_back(x);
        ys.push_back(y);
        zs.push_back(z);

        current = model->next_monomers[current];

    }

    model->count_X[abs(xs.back())]+=1;
    model->count_Y[abs(ys.back())]+=1;

    long double r_g = 0;
    long double xdiff1, ydiff1, zdiff1;
    double  A_element=0, D_element=0, BC_element2 = 0;

    //long long int second_current;

    int number_of_monomers = model->number_of_spins();

    for (int e = 0; e < number_of_monomers; e++)
    {
        for (int e1 = 0; e1 < number_of_monomers; e1++) {
            //расстояние на торе
            xdiff1 = xs[e1] - xs[e];
            ydiff1 = ys[e1] - ys[e];
            zdiff1 = zs[e1] - zs[e];

            r_g = r_g + xdiff1 * xdiff1 + ydiff1 * ydiff1 + zdiff1*zdiff1;

            A_element +=xdiff1 * xdiff1;
            D_element += ydiff1 * ydiff1;
            BC_element2 += xdiff1*ydiff1;

        }

    }

    model->gyration << 0.5*r_g/number_of_monomers/number_of_monomers;
    A_element = 1.0*A_element/number_of_monomers/number_of_monomers/2.0;
    D_element = 1.0*D_element/number_of_monomers/number_of_monomers/2.0;
    BC_element2 = 1.0*BC_element2/number_of_monomers/number_of_monomers/2.0;

    double D = (A_element+D_element)*(A_element+D_element)-4*(A_element*D_element - BC_element2*BC_element2);
    double eig1 = ((A_element+D_element) + sqrt(D))*0.5;
    double eig2 =((A_element+D_element) - sqrt(D))*0.5;

    model->eigval1 << eig1;
    model->eigval2 << eig2;
    model->aratio << 1.0*(eig1-eig2)*(eig1-eig2)/((eig1+eig2)*(eig1+eig2));

}



void Lattice_2D::radius(SAW_model *model) {
    long int point1x = model->end_conformation % lattice_side;
    long int point1y = model->end_conformation / lattice_side;
    long int point1xs = model->start_conformation % lattice_side;
    long int point1ys = model->start_conformation / lattice_side;

    //расстояние на торе
    long int xdiff = abs(point1x- point1xs);
    if (xdiff > (lattice_side  / 2))
        xdiff = lattice_side - xdiff;

    long int ydiff = abs(point1y- point1ys);
    if (ydiff > (lattice_side / 2))
        ydiff = lattice_side - ydiff;


    long int r = xdiff *xdiff  + ydiff*ydiff ;
    model->dists << r;

    model->count_R2[r]+=1;
    //model->count_X[xdiff]+=1;
    //model->count_Y[ydiff]+=1;

    std::vector <long long int> xs,ys;
    xs.push_back(0);
    ys.push_back(0);

    int direction;


    long int current = model->start_conformation;

    //std::vector <std::vector<int>> steps = {   {1,0}, {-1,0}, {0,1}, {0,-1} };
    for (long long int i =1; i < model->number_of_spins(); i++)
    {
        direction = model->directions[current];
        long long int x = xs.back()+steps[direction][0];
        long long int y = ys.back()+steps[direction][1];

        xs.push_back(x);
        ys.push_back(y);


        current = model->next_monomers[current];

    }

    model->count_X[abs(xs.back())]+=1;
    model->count_Y[abs(ys.back())]+=1;

    long double r_g = 0;
    long double xdiff1, ydiff1;
    double  A_element=0, D_element=0, BC_element2 = 0;

    //long long int second_current;

    int number_of_monomers = model->number_of_spins();

    for (int e = 0; e < number_of_monomers; e++)
    {
        for (int e1 = 0; e1 < number_of_monomers; e1++) {
            //расстояние на торе
            xdiff1 = xs[e1] - xs[e];
            ydiff1 = ys[e1] - ys[e];

            r_g = r_g + xdiff1 * xdiff1 + ydiff1 * ydiff1;

            A_element +=xdiff1 * xdiff1;
            D_element += ydiff1 * ydiff1;
            BC_element2 += xdiff1*ydiff1;

        }

    }

    model->gyration << 0.5*r_g/number_of_monomers/number_of_monomers;
    A_element = 1.0*A_element/number_of_monomers/number_of_monomers/2.0;
    D_element = 1.0*D_element/number_of_monomers/number_of_monomers/2.0;
    BC_element2 = 1.0*BC_element2/number_of_monomers/number_of_monomers/2.0;

    double D = (A_element+D_element)*(A_element+D_element)-4*(A_element*D_element - BC_element2*BC_element2);
    double eig1 = ((A_element+D_element) + sqrt(D))*0.5;
    double eig2 =((A_element+D_element) - sqrt(D))*0.5;

    model->eigval1 << eig1;
    model->eigval2 << eig2;
    model->aratio << 1.0*(eig1-eig2)*(eig1-eig2)/((eig1+eig2)*(eig1+eig2));

}

