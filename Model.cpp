//
// Created by kamilla on 01.02.2022.
//
//#include "Lattice.h"

#include <iostream>
#include <fstream>

#include <queue>
#include "Model.h"


std::random_device generator;
std::uniform_real_distribution<double> distribution(0.0,1.0);
//std::mt19937 generator(1234);
std::random_device generators1;

std::uniform_int_distribution<int> distribution2(0, 1);
//std::mt19937 generators3(12377);
std::random_device generators3;

std::uniform_real_distribution<double> distribution_theta(-PI, PI);
std::random_device generatorstheta;

//change after
std::uniform_int_distribution<int> distribution1(0, 3);
//std::mt19937 generators1(123);
std::random_device generators2;


SAW_model::SAW_model(short n, short d) {
    L = n;
    dimension = d;
    //create lattice
    if (d == 2)
    {
        lattice = new Lattice_2D(n+5);
    }
    else if (d == 3)
    {
        lattice = new Lattice_3D(n+5);
    }
    else
    {
        return ;
    }

    if (lattice== nullptr) return ;

    // 0 - отстутствие элемента на решетке
    sequence_on_lattice.resize(lattice->NumberOfNodes(),-5); //последовательность мономеров
    next_monomers.resize(lattice->NumberOfNodes(),-1 ); //номер-ссылка на следующий узел
    previous_monomers.resize(lattice->NumberOfNodes(), -1 ); //номер-ссылка на предыдующий узел
    directions.resize(lattice->NumberOfNodes(),-1 ); //направления из {0,1,2,3}

}

void SAW_model::save_geometry() {
    lattice->radius(this);
}

void SAW_model::calc_bulc() {
    double bulk4_now=0, bulk3_now=0,bulk2_now=0,bulk5_now=0, bulk6_now=0;;

    int number_of_monomers = L;
    long int current = start_conformation;
    long int step;
    int k = 0;
    for (int e = 0; e < number_of_monomers; e++)
    {
        k = 0;
        for (int j = 0; j < lattice->ndim2(); j++) {
            step = lattice->map_of_contacts_int[lattice->ndim2() * current + j];
            if (sequence_on_lattice[step] != -5) {
                k +=1;
            }
        }

        if(k==2) {
            bulk2_now+=1;
        }
        if(k==3) {
            bulk3_now+=1;
        }
        if(k==4) {
            bulk4_now+=1;
        }
        if (dimension==3) {
            if (k == 5) {
                bulk5_now += 1;
            }
            if (k == 6) {
                bulk6_now += 1;
            }
        }

        current = next_monomers[current];
    }
    bulk2 << 1.0*bulk2_now/number_of_monomers;
    bulk3 << 1.0*bulk3_now/number_of_monomers;
    bulk4 << 1.0*bulk4_now/number_of_monomers;
    if (dimension==3) {
        bulk5 << 1.0*bulk5_now/number_of_monomers;
        bulk6 << 1.0*bulk6_now/number_of_monomers;
    }

}



void XY_SAW::FlipMove(double J)
{


    long double new_E = 0.;
    double hh = 0.;
    short rand_path = distribution2(generators3);
    int step_on_lattice;
    long int step;
    coord_t new_point;
    double oldspin;
    long int del, temp;

    if (rand_path == 0) {//переставляем начало в конец

        step_on_lattice = distribution1(generators1);
        new_point = lattice->map_of_contacts_int[lattice->ndim2() * end_conformation + step_on_lattice];
        oldspin = sequence_on_lattice[start_conformation];

        if (sequence_on_lattice[new_point] == -5) { //проверка, что в узле нет мономеров

            //делаем апдейт

            //удаляем начало
            temp = start_conformation;
            start_conformation = next_monomers[start_conformation];
            next_monomers[temp] = -1;
            previous_monomers[start_conformation] = -1;
            sequence_on_lattice[temp] = -5;
            //смотрим потери
            for (int j = 0; j < lattice->ndim2(); j++) {
                step = lattice->map_of_contacts_int[lattice->ndim2() * temp + j];
                if (sequence_on_lattice[step] != -5) {
                    hh = hh - cos(oldspin - sequence_on_lattice[step]);
                }
            }

            //добавляем в конец
            next_monomers[end_conformation] = new_point;
            sequence_on_lattice[new_point] =  distribution_theta(generatorstheta); //выбор спина
            previous_monomers[new_point] = end_conformation;
            end_conformation = new_point;
            //смотрим выигрыш
            for (int j = 0; j < lattice->ndim2(); j++) {
                step = lattice->map_of_contacts_int[lattice->ndim2() * end_conformation + j];
                if (sequence_on_lattice[step] != -5) {
                    hh = hh + cos(sequence_on_lattice[end_conformation] - sequence_on_lattice[step]);
                }
            }

            new_E = E + hh;
            //new_H = current_H_counts + sequence_on_lattice[new_point] - sequence_on_lattice[start_conformation];
            //new_H = current_H_counts + sequence_on_lattice[new_point] - sequence_on_lattice[temp];

            //p1 = exp( -(-(new_E - E) * J - (new_H - current_H_counts) * h));
            double p1 = exp( -( -(new_E - E) * J ));
            double p_metropolis = std::min(1.0, p1);
            double q_rd = distribution(generator);
            if (q_rd < p_metropolis) { //принимаем изменения
                E = new_E;
                //current_H_counts = new_H;
                sequence_on_lattice[temp] = -5; //делаю здесь, так как проще считать энергию(!!!)

                //корректируем информацию о направлениях
                directions[temp] = -1;
                directions[previous_monomers[end_conformation]] = step_on_lattice;

            }
            else {//отменяем изменения
                //удаляем конец
                del = end_conformation;
                end_conformation = previous_monomers[end_conformation];
                next_monomers[end_conformation] = -1;
                previous_monomers[del] = -1;
                sequence_on_lattice[del] = -5;

                //возвращаем начало
                previous_monomers[start_conformation] = temp;
                next_monomers[temp] = start_conformation;
                start_conformation = temp;
                sequence_on_lattice[start_conformation] = oldspin;

            }

        }
        else
        {
            //места нет, выходим из шага

        }

    }
    else {//переставляем конец в начало

        step_on_lattice = distribution1(generators1);
        new_point = lattice->map_of_contacts_int[lattice->ndim2() * start_conformation + step_on_lattice];
        oldspin = sequence_on_lattice[end_conformation];

        if (sequence_on_lattice[new_point] == -5) { //проверка, что в узле нет мономеров

            //делаем апдейт

            //удаляем конец
            temp = end_conformation;
            end_conformation = previous_monomers[end_conformation];
            /*if (previous_monomers[end_conformation] < 0) {
                std::cout << "problem update " << std::endl;
            }*/
            previous_monomers[temp] = -1;
            next_monomers[end_conformation] = -1;
            sequence_on_lattice[temp] = -5;
            //смотрим потери
            for (int j = 0; j < lattice->ndim2(); j++) {
                step = lattice->map_of_contacts_int[lattice->ndim2() * temp + j];
                if (sequence_on_lattice[step] != -5) {
                    hh = hh - cos(oldspin -sequence_on_lattice[step]);
                }
            }

            //добавляем в начало
            previous_monomers[start_conformation] = new_point;
            sequence_on_lattice[new_point] = distribution_theta(generatorstheta); //выбор спина
            next_monomers[new_point] = start_conformation;
            start_conformation = new_point;
            //смотрим выигрыш
            for (int j = 0; j < lattice->ndim2(); j++) {
                step = lattice->map_of_contacts_int[lattice->ndim2() * start_conformation + j];
                if (sequence_on_lattice[step] != -5) {
                    hh = hh + cos(sequence_on_lattice[start_conformation] - sequence_on_lattice[step]);
                }
            }

            new_E = E + hh;
            //new_H = current_H_counts + sequence_on_lattice[new_point] - sequence_on_lattice[temp];

            //p1 = exp(-(new_E - E) * J - (new_H - current_H_counts) * h);

            double p1 = exp( -( -(new_E - E) * J ));
            double p_metropolis = std::min(1.0, p1);
            double q_rd = distribution(generator);

            if (q_rd < p_metropolis) {
                E = new_E;
                //current_H_counts = new_H;
                sequence_on_lattice[temp] = -5; //делаю здесь, так как проще считать энергию(!!!)

                //корректируем информацию о направлениях
                directions[end_conformation] = -1;
                directions[start_conformation] = lattice->inverse_steps[step_on_lattice];

            }
            else {//отменяем изменения
                //удаляем начало
                del = start_conformation;
                start_conformation = next_monomers[start_conformation];
                previous_monomers[start_conformation] = -1;
                next_monomers[del] = -1;
                sequence_on_lattice[del] = -5;

                //возвращаем конец
                next_monomers[end_conformation] = temp;
                previous_monomers[temp] = end_conformation;
                end_conformation = temp;
                sequence_on_lattice[end_conformation] = oldspin;


            }
        }
        else {
            //некуда идти
        }

    }

}

void XY_SAW ::ClusterStep(double J)
{
    // делаем кластерный апдейт
    std::uniform_int_distribution<long int> distribution_spin(0, L-1);
    std::random_device generator_spin;
    long int choose_spin = distribution_spin(generator_spin);

    long int coord = start_conformation;
    for (long int spin = 1; spin < choose_spin; spin++)
    {
        coord = next_monomers[coord];
    }

    double flipdirection = distribution_theta(generatorstheta);

    double x1 = cos(sequence_on_lattice[coord])*cos(flipdirection)+sin(sequence_on_lattice[coord])*sin(flipdirection);
    double s = sin(sequence_on_lattice[coord])-2*x1*sin(flipdirection);
    double c = cos(sequence_on_lattice[coord])-2*x1*cos(flipdirection);

    if (s<-1.) s=-1;
    if (s>1.) s=1;
    if (c<-1.) c=-1;
    if (c>1.) c=1;

    sequence_on_lattice[coord] = (s > 0) ? acos(c) : -acos(c);
    int sign = (x1 < 0) ? -1 : (x1 > 0);
    double x = x1;
    //std::valarray<bool> used_coords;
    //used_coords.resize(lattice->NumberOfNodes(), false  );

    std::queue<long int> Cluster;

    Cluster.push(coord);
    used_coords[coord] = true;
    double tempscalar = 0;
    int tempsign = -1;
    long int temp, step;
    while (!Cluster.empty()) {
        temp = Cluster.front();
        Cluster.pop();

        for (int j = 0; j < lattice->ndim2(); j++)
        {
            step = lattice->map_of_contacts_int[lattice->ndim2() * temp + j];
            tempscalar = cos(sequence_on_lattice[step])*cos(flipdirection)+sin(sequence_on_lattice[step])*sin(flipdirection);
            tempsign = (tempscalar < 0) ? -1 : (tempscalar > 0);

            double P_add =  1 - exp(-2*J*tempscalar*x);

            double p = distribution(generator);
            //???
            if ( sequence_on_lattice[step]!=-5. &&
                 tempsign == sign &&
                 p < P_add &&
                 !used_coords[step]) {
                Cluster.push(step);
                used_coords[step]= true;

                double s = sin(sequence_on_lattice[step])-2*tempscalar*sin(flipdirection);
                double c = cos(sequence_on_lattice[step])-2*tempscalar*cos(flipdirection);
                if (s<-1.) s=-1;
                if (s>1.) s=1;
                if (c<-1.) c=-1;
                if (c>1.) c=1;
                sequence_on_lattice[step] = (s > 0) ? acos(c) : -acos(c);
            }
        }
    }

    long int current = start_conformation;
    for (int e = 0; e < L; e++)
    {
        used_coords[current] = false;

        current = next_monomers[current];
    }

    Energy();
}

void SAW_model::Reconnect(int j)
{
    long int c;
    long int step = lattice->map_of_contacts_int[lattice->ndim2()*end_conformation +j];

    //проверка, что проверенный узел занят спином
    if (sequence_on_lattice[step] == -5 || next_monomers[step] == -1 ||
            step == previous_monomers[end_conformation]) {
        return;
    }

    long int new_end = next_monomers[step];

    next_monomers[step]=end_conformation;
    directions[step]=lattice->inverse_steps[j];
    c = end_conformation;
    long int new_c;
    while (c!=new_end)
    {
        new_c=previous_monomers[c];
        next_monomers[c]=previous_monomers[c];
        directions[c]=lattice->inverse_steps[directions[new_c]];
        c=new_c;
    }
    long int temp_prev_next = next_monomers[new_end];
    previous_monomers[end_conformation]=step;
    c=end_conformation;
    while (c!=new_end)
    {
        new_c=next_monomers[c];
        previous_monomers[new_c]=c;
        c=new_c;
    }
    end_conformation=new_end;
    previous_monomers[new_end]=temp_prev_next;
    next_monomers[new_end]=-1;
    directions[new_end]=-1;
}

XY_SAW::XY_SAW(short int n, short d) : SAW_model(n, d) {
    start_conformation=0;
    end_conformation=n-1;

    for (int i = 1; i < n-1; i++)
    {
        previous_monomers[i]=i-1;
        sequence_on_lattice[i]=PI;
        next_monomers[i]=i+1;
    }
    sequence_on_lattice[0] = PI;
    sequence_on_lattice[end_conformation] = PI; //начальная последовательность
    next_monomers[0] = 1;
    previous_monomers[n-1] = n-2;
    E =  -(n-1);

    //сначала все направления - движение вправо
    for (int i = 0; i < n-1; i++)
    {
        directions[i]=0;
    }

    used_coords.resize(lattice->NumberOfNodes(), false  );
}

void XY_SAW::Energy() {
    double hh = 0;
    long int current_position = start_conformation;
    coord_t  step;
    long int mag = 0;
    for (int i =0; i<L; i++){
        for ( int j=0; j<lattice->ndim2(); j++ ){
            step = lattice->map_of_contacts_int[lattice->ndim2()*current_position+j];
            if ( sequence_on_lattice[step]!=-5  )
            {
                hh=hh+cos(sequence_on_lattice[current_position]-sequence_on_lattice[step]);
            }
        }
        current_position=next_monomers[current_position];
    }
    E = -(hh/2.0);
}

void XY_SAW::save_measurements() {

    int number_of_monomers = L;

    energy << 1.0*(E)/number_of_monomers;
    energy_sq << 1.0*(E)/number_of_monomers* 1.0*(E)/number_of_monomers;
    energy_4 << 1.0*(E)/number_of_monomers* 1.0*(E)/number_of_monomers* 1.0*(E)/number_of_monomers* 1.0*(E)/number_of_monomers;

    long int current = start_conformation;
    double sum_sin_2 = 0.0;
    double sum_cos_2 = 0.0;
    double sum_sin_1 = 0.0;
    double sum_cos_1 = 0.0;
    for (int e = 0; e < number_of_monomers; e++)
    {
        sum_sin_1  += sin(sequence_on_lattice[current]);
        sum_cos_1  += cos(sequence_on_lattice[current]);
        sum_sin_2  += (sin(sequence_on_lattice[current])*sin(sequence_on_lattice[current]) );
        sum_cos_2  += (cos(sequence_on_lattice[current])*cos(sequence_on_lattice[current]) );
        current = next_monomers[current];
    }

    sum_sin_1/=number_of_monomers;
    sum_cos_1/=number_of_monomers;
    sum_sin_2/=number_of_monomers;
    sum_cos_2/=number_of_monomers;
    sum_sin_2/=number_of_monomers;
    sum_cos_2/=number_of_monomers;

    mags_sin << sum_sin_1;
    mags_cos << sum_cos_1;

    magnetization_sq <<sum_sin_1*sum_sin_1 + sum_cos_1*sum_cos_1;
    magnetization_4 << (sum_sin_1*sum_sin_1 + sum_cos_1*sum_cos_1)*(sum_sin_1*sum_sin_1 + sum_cos_1*sum_cos_1);

    int N_cos = (sum_cos_1-(-1.))/h_l;
    count_cos[N_cos]+=1;
    int N_m2 = ((sum_sin_1*sum_sin_1 + sum_cos_1*sum_cos_1)-(-1.))/h_l;
    count_m2[N_m2]+=1;
    int N_E = ((1.0*(E)/number_of_monomers)-(-2.))/h_l;
    count_E[N_E]+=1;
}

void XY_SAW::write_file(long int i, double J) {
    std::string filename;
    std::ofstream out_result;


    //!!!!!!!!!!!!!!!!!! temporary
    double h=0;

    filename = "XY_"+std::to_string(J)+"_"+std::to_string(h)+"_"+std::to_string(L)+"_"+std::to_string(0)+".txt";

    out_result.open(filename);
    //out_result << mc_steps<<" " << number_of_monomers << " " << J << " " << h  <<   " ";
    out_result << "N J h mean_R_sq err_mean_R_sq mean_R_gyr_sq err_mean_R_gyr_sq ";
    out_result << "mean_e err_mean_e mean_e_sq err_mean_e_sq mean_e_fourth err_mean_e_fourth ";

    out_result << "mean_sin err_mean_sin mean_cos err_mean_cos mean_m_sq err_mean_m_sq mean_m_fourth err_mean_m_fourth steps " ;
    if (dimension==3) out_result << "bulk2 err_bulk2 bulk3 err_bulk3 bulk4 err_bulk4 bulk5 err_bulk5 bulk6 err_bulk6 ";
    else out_result << "bulk2 err_bulk2 bulk3 err_bulk3 bulk4 err_bulk4 ";


    out_result << "lambda1 err_lambda1 lambda2 err_lambda2 asperical err_aspherical "<< std::endl;

    out_result << L << " " << J << " " << h <<  " ";
    out_result << dists.mean() << " " << dists.errorbar()<< " " << gyration.mean() << " " << gyration.errorbar() << " ";

    out_result << energy.mean() << " " << energy.errorbar() << " ";
    out_result << energy_sq.mean() << " " << energy_sq.errorbar() << " ";
    out_result << energy_4.mean() << " " << energy_4.errorbar() << " ";

    out_result << mags_sin.mean() << " " << mags_sin.errorbar() << " ";
    out_result << mags_cos.mean() << " " << mags_cos.errorbar() << " ";
    out_result << magnetization_sq.mean() << " " << magnetization_sq.errorbar() << " ";
    out_result << magnetization_4.mean() << " " << magnetization_4.errorbar() << " ";
    out_result << i << " ";

    out_result << bulk2.mean() << " " << bulk2.errorbar() << " ";
    out_result << bulk3.mean() << " " << bulk3.errorbar() << " ";
    out_result << bulk4.mean() << " " << bulk4.errorbar() << " ";
    if (dimension==3) {
        out_result << bulk5.mean() << " " << bulk5.errorbar() << " ";
        out_result << bulk6.mean() << " " << bulk6.errorbar() << " ";
    }

    out_result << eigval1.mean() << " " << eigval1.errorbar() << " ";
    out_result << eigval2.mean() << " " << eigval2.errorbar() << " ";
    out_result << aratio.mean() << " " << aratio.errorbar() << " ";


    out_result << std::endl;

    out_result.close();


    out_result.close();


    filename = "R2_"+std::to_string(J)+"_"+std::to_string(h)+"_"+std::to_string(L)+".txt";

    //filename = "Radius_"+std::to_string(J)+"_"+std::to_string(number_of_monomers)+"_CanonicalIsing.txt";

    out_result.open(filename);
    //out_result << mc_steps<<" " << number_of_monomers << " " << J << " " << h  <<   " ";

    out_result << "N J h mean_R_sq err_mean_R_sq mean_R_gyr_sq err_mean_R_gyr_sq " << std::endl;

    out_result << L << " " << J << " " << h <<  " ";
    out_result << dists.mean() << " " << dists.errorbar()<< " " << gyration.mean() << " " << gyration.errorbar() << std::endl;

    for (auto c : count_R2)
    {
        out_result << c.first << " " << c.second << std::endl;
    }

    out_result.close();

    filename = "X_"+std::to_string(J)+"_"+std::to_string(h)+"_"+std::to_string(L)+".txt";


    //filename = "Radius_"+std::to_string(J)+"_"+std::to_string(number_of_monomers)+"_CanonicalIsing.txt";

    out_result.open(filename);
    //out_result << mc_steps<<" " << number_of_monomers << " " << J << " " << h  <<   " ";

    out_result << "N J h mean_R_sq err_mean_R_sq mean_R_gyr_sq err_mean_R_gyr_sq " << std::endl;

    out_result << L << " " << J << " " << h <<  " ";
    out_result << dists.mean() << " " << dists.errorbar()<< " " << gyration.mean() << " " << gyration.errorbar() << std::endl;

    for (auto c : count_X)
    {
        out_result << c.first << " " << c.second << std::endl;
    }

    out_result.close();

    filename = "Y_"+std::to_string(J)+"_"+std::to_string(h)+"_"+std::to_string(L)+".txt";


    //filename = "Radius_"+std::to_string(J)+"_"+std::to_string(number_of_monomers)+"_CanonicalIsing.txt";

    out_result.open(filename);
    //out_result << mc_steps<<" " << number_of_monomers << " " << J << " " << h  <<   " ";

    out_result << "N J h mean_R_sq err_mean_R_sq mean_R_gyr_sq err_mean_R_gyr_sq " << std::endl;

    out_result << L << " " << J << " " << h <<  " ";
    out_result << dists.mean() << " " << dists.errorbar()<< " " << gyration.mean() << " " << gyration.errorbar() << std::endl;

    for (auto c : count_Y)
    {
        out_result << c.first << " " << c.second << std::endl;
    }

    out_result.close();


    filename = "Counts_E_"+std::to_string(J)+"_"+std::to_string(h)+"_"+std::to_string(L)+".txt";


    //filename = "Radius_"+std::to_string(J)+"_"+std::to_string(number_of_monomers)+"_CanonicalIsing.txt";

    out_result.open(filename);
    //out_result << mc_steps<<" " << number_of_monomers << " " << J << " " << h  <<   " ";

    out_result << "N J h mean_R_sq err_mean_R_sq mean_R_gyr_sq err_mean_R_gyr_sq " << std::endl;

    out_result << L << " " << J << " " << h <<  " ";
    out_result << dists.mean() << " " << dists.errorbar()<< " " << gyration.mean() << " " << gyration.errorbar() << std::endl;

    for (auto c : count_E)
    {
        out_result << c.first << " " << c.second << std::endl;
    }

    out_result.close();

    filename = "Counts_cos_"+std::to_string(J)+"_"+std::to_string(h)+"_"+std::to_string(L)+".txt";


    //filename = "Radius_"+std::to_string(J)+"_"+std::to_string(number_of_monomers)+"_CanonicalIsing.txt";

    out_result.open(filename);
    //out_result << mc_steps<<" " << number_of_monomers << " " << J << " " << h  <<   " ";

    out_result << "N J h mean_R_sq err_mean_R_sq mean_R_gyr_sq err_mean_R_gyr_sq " << std::endl;

    out_result << L << " " << J << " " << h <<  " ";
    out_result << dists.mean() << " " << dists.errorbar()<< " " << gyration.mean() << " " << gyration.errorbar() << std::endl;

    for (auto c : count_cos)
    {
        out_result << c.first << " " << c.second << std::endl;
    }

    out_result.close();



    filename = "Counts_mag2_"+std::to_string(J)+"_"+std::to_string(h)+"_"+std::to_string(L)+".txt";


    //filename = "Radius_"+std::to_string(J)+"_"+std::to_string(number_of_monomers)+"_CanonicalIsing.txt";

    out_result.open(filename);
    //out_result << mc_steps<<" " << number_of_monomers << " " << J << " " << h  <<   " ";

    out_result << "N J h mean_R_sq err_mean_R_sq mean_R_gyr_sq err_mean_R_gyr_sq " << std::endl;

    out_result << L << " " << J << " " << h <<  " ";
    out_result << dists.mean() << " " << dists.errorbar()<< " " << gyration.mean() << " " << gyration.errorbar() << std::endl;

    for (auto c : count_m2)
    {
        out_result << c.first << " " << c.second << std::endl;
    }

    out_result.close();


}




