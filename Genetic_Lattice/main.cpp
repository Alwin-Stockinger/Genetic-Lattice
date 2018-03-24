#include<iostream>
#include<random>
#include<vector>
#include<cmath>
#include<iostream>
#include<algorithm>

#define PI 3.14159265


using namespace std;

int const a=10;
int n=100;

int kmax=20;
int lmax=20;

int parent_amount=2;        //changes in code have to be made!!!
double mutation_rate=0.01;


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
                    sum+=v(sqrt(pow(k+l*x*cos(phi),2)+pow(l*x*sin(phi),2)));            //optimizable
                }
            }
        }
        return sum;
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
            //cout<<sum<<endl<<tric_sum<<endl;
            return exp(1-sum/tric_sum);
        }
      
        

} fit;


class Individual{
    vector<bool> x;
    vector<bool> phi;
    double fitness;
    

    public:

    vector<bool> getVecX(){
        return x;
    }
    vector<bool> getVecPhi(){
        return phi;
    }


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
        calcFitness();
        return fitness;
    }

    void calcFitness(){
        fitness=fit(getX(),getPhi());
    }

    //Inital Creator
    Individual(int x_length,int phi_length){
        random_device rd;
        mt19937 gen(rd());
        uniform_int_distribution<> dis(0,1);
    
        for(int i=0;i<x_length;i++) x.push_back(dis(gen));
        for(int i=0;i<phi_length;i++) phi.push_back(dis(gen));
    }
    //child Creator
    Individual(vector<bool> xVec, vector<bool> phiVec){ //has to be changed for more or less than 2 parents!!!!!!
        x=xVec;
        phi=phiVec;
        mutate();
    }

    private:
    void mutate(){
        random_device rd;
        mt19937 gen(rd());
        uniform_real_distribution<> dis(0,1);
        for(int i=0;i<x.size();++i){
            if(dis(gen)<=mutation_rate){
                x.at(i)=!x.at(i);
            }
        }
        for(int i=0;i<phi.size();++i){
            if(dis(gen)<=mutation_rate){
                phi.at(i)=!phi.at(i);
            }
        }
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
        //cout<<"children will now be born"<<endl;
        for(int i=0;i<size;i+=2){      // double because 2 children are created
            vector<Individual> child=genChild();
            cout<<"First child generated"<<endl;
            children.insert(children.end(),child.begin(),child.end());
        }
        //cout<<"Children become Parents"<<endl;
        indiv=children;
    }

    void printIndiv(){
        auto print=[](Individual ind){cout<<"Individual:"<<endl<<"X="<<ind.getX()<<endl<<"Phi="<<ind.getPhi()<<endl<<"Fit="<<ind.getFitness()<<endl<<endl;};
        for_each(indiv.begin(),indiv.end(),print);
    }


    private:

    Individual getParent(){
        random_device rd;
        mt19937 gen(rd());
        uniform_real_distribution<> dis(0,genFitness);     
        double rand=dis(gen);
        int i=0;
        while(indiv.at(i).getFitness()<rand){
            rand-=indiv.at(i).getFitness();
            i++;
        }
        return indiv.at(i);
    }


    vector<Individual> getParents(int amount){
        vector<Individual> parents;
        for(int i=0;i<amount;++i) parents.push_back(getParent());
        return parents;
    }

    vector<Individual> genChild(){  //has to be changed for more or less than two parents!
        vector<Individual> parents=getParents(parent_amount);
        //cout<<"Parents are now known!"<<endl;
        random_device rd;
        mt19937 gen(rd());
        uniform_int_distribution<> disX(0,parents.at(0).getVecX().size());
        uniform_int_distribution<> disPhi(0,parents.at(0).getVecPhi().size());

        int cutX=disX(gen);
        int cutPhi=disPhi(gen);

        vector<bool> x1(parents.at(0).getVecX().size());
        vector<bool>phi1(parents.at(0).getVecPhi().size());
        vector<bool>x2(parents.at(0).getVecX().size());
        vector<bool>phi2(parents.at(0).getVecPhi().size());

        //cout<<"Genes will be created!"<<endl;
        for(int i=0;i<cutX;++i){
            x1[i]=parents[0].getVecX().at(i);
            x2[i]=parents[1].getVecX().at(i);
        }
        //cout<<"Secon Gene now edited"<<endl;
        for(int i=cutX;i<parents.at(0).getVecX().size();++i){
            x1[i]=parents[1].getVecX().at(i);
            x2[i]=parents[0].getVecX().at(i);
        }

        for(int i=0;i<cutPhi;++i){
            phi1[i]=parents[0].getVecPhi().at(i);
            phi2[i]=parents[1].getVecPhi().at(i);
        }
        for(int i=cutPhi;i<parents.at(0).getVecPhi().size();++i){
            phi1[i]=parents[1].getVecPhi().at(i);
            phi2[i]=parents[0].getVecPhi().at(i);
        }

        //cout<<"The two children will be created"<<endl;
        Individual child1(x1,phi1);
        Individual child2(x2,phi2);

        vector<Individual> children;
        children.push_back(child1);
        children.push_back(child2);
        return children;
    }

    void calcGenFitness(){
        genFitness=accumulate(indiv.begin(),indiv.end(),0.,
            [](double x, Individual y){ 
            return x+y.getFitness();
            });
    }

};




int main(){



    Generation gen(8,8,8);

    gen.printIndiv();
    cout<<endl<<"Next Generation"<<endl;
    gen.nextGen(8);
    
    gen.printIndiv();
    bool x;
    cin>>x;
    return 0;
}