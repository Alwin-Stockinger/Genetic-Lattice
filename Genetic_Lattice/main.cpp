#include<iostream>
#include<random>
#include<vector>
#include<cmath>
#include<iostream>

#define PI 3.14159265


using namespace std;

int const a=10;
int n=100;

int kmax=20;
int lmax=20;



struct Tric{
    public:
    double const phi=PI/3.;
    double const x=1.;
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
        
        double operator() (double x,double phi) const{
            double sum=latticeSum(x,phi,potential);
            cout<<sum<<endl<<tric_sum<<endl;
            return exp(1-sum/tric_sum);
        }
      
        

};


class Individual{
    vector<bool> x;
    vector<bool> phi;
    double fitness;
    Fitness fit;

    public:
    double getX(){
        double accum=accumulate(x.rbegin(),x.rend(),0,[](int i, int j){ return (i<<1)+j;});    //shift the x bits to the left, so that y can be inserted on right site
        
        accum/=pow(2,x.size());
        
        return accum;
    }
    double getPhi(){
        double accum=accumulate(phi.rbegin(),phi.rend(),0,[](int x, int y){ return (x<<1)+y;});
        accum/=pow(2,phi.size())*PI/2.;
        return accum;
    }
    double getFitness(){
        return fitness;
    }

    void calcFitness(){
        fitness=fit(getX(),getPhi());
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

    double genFitness=0;

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

    void nextGen(int size){
        calcGenFitness();

        vector<Individual> children;
        for(int i=0;i<size;++i){
            children.push_back(genChild());
        }
    }

    private:
    Individual genChild(){
        
    }

    void calcGenFitness(){
        double sum=accumulate(indiv.begin(),indiv.end(),0.,
            [](double x, Individual y){ 
            return x+y.getFitness();
            });
        genFitness=sum;
    }

};




int main(){
    Fitness fitness;
    Individual ind(7,7);
    ind.calcFitness();

    cout<<ind.getFitness()<<endl;
    int x;
    cin>>x;
    return 0;
}