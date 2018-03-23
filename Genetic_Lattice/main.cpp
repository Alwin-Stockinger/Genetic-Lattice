#include<iostream>
#include<random>
#include<vector>
#include<cmath>


using namespace std;

int const a=1;
int n=100;

int kmax=20;
int lmax=20;

class Individual;


struct Tric{
    public:
    double const phi=60;
    double const x=a;
} tric;

class Potential{
    public:
    double operator() (double r) const{
        return exp(-r)/r;
    }
};

class LatticeSum{
    public:
    double operator() (double x, double phi, Potential v) const{
        double sum=0;
        for(int k=-kmax;k<=kmax;++k){
            for(int l=-lmax;l<=lmax;++l){
                if(l!=0||k!=0){
                    sum+=v(sqrt(k*k+2*k*l*sin(phi)+l*l*x*x));
                }
            }
        }
    }
};

class Fitness{

    
    LatticeSum latticeSum;
    Potential potential;
    public:
    double tric_sum;


        Fitness(){
            tric_sum=latticeSum(tric.x,tric.phi,potential);
        }
    
        double operator() (Individual ind) const{
            return 
        }
    
};


class Individual{
    vector<bool> x;
    vector<bool> phi;
    double fitness;

    public:
    vector<bool> getX(){
        return x;
    }
    vector<bool> getPhi(){
        return phi;
    }


    Individual(int x_length,int phi_length){
        random_device rd;
        mt19937 gen(rd());
        uniform_int_distribution<> dis(0,1);
    
        for(int i=0;i<x_length;i++) x.push_back(dis(gen));
        for(int i=0;i<phi_length;i++) phi.push_back(dis(gen));
    }


    
    
};

class Generation{
    vector<Individual> indiv;
    public:
    vector<Individual> getIndividuals(){
        return indiv;
    } 

    Generation(int ind_amount, int x_length, int phi_length){
        for(int i=0;i<ind_amount;++i){
            Individual ind(x_length,phi_length);
            indiv.push_back(ind);
        }
    }

    private:
  

};




int main(){
    Fitness fit;
    cout<<fit.tric_sum;
    int x;
    cin>>x;
    return 0;
}