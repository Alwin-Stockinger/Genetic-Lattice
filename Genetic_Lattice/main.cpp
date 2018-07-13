#include<iostream>
#include<random>
#include<vector>
#include<cmath>
#include<iostream>
#include<algorithm>
#include<string>
#include<bitset>

#include<thread>
#include<mutex>
#include<chrono>
#include<ctime>

#include "matplotlibcpp.h"

#define PI 3.141592653589793238462643383279502884


using namespace std;


double const lambda=1;
double d=lambda*10;
int n=100;
double gen_i=1;
int const bitlen=10;


int layers=2;
double volDensity=0.1;      
double density=d*volDensity/layers;  //density in Layer



double const rcut=lambda*6;
//int kmax=15;
//int lmax=15;

int parent_amount=2;        //changes in code have to be made!!!
double mutation_rate=0.05;

vector<double> veczero={0,0,0};

int threadcount=4;

int survivors=1;//this many best genes will deifnitly survive

struct Tric{
    public:
    double const phi=PI/2.;
    double const x=1;
    double a;
    vector<double> cx;
    vector<double> cy;
    vector<double> hTric;

    Tric(){
        a=1./sqrt(x*sin(phi)*density);
        setH();
        setCx();
        setCy();
    }
    void calculateTric(){
        a=1./sqrt(x*sin(phi)*density);
        setH();
        setCx();
        setCy();
    }

    void setH(){
        hTric.clear();
        for(int i=0;i<layers-1;i++){
            hTric.push_back(d/(layers-1.));
        }
    }

    void setCx(){
        cx.clear();
        for(int i=0;i<layers-1;i++){
            if(i%2) cx.push_back(a/(2.));
            else cx.push_back(a/(-2.));
        }
    }

    void setCy(){
        cy.clear();
        for(int i=0;i<layers-1;i++){
            if(i%2) cy.push_back(a*x*sin(phi)/(2.));
            else cy.push_back(-a*x*sin(phi)/(2.));
        }
    }



} tric;

class Potential{
    public:
    double operator() (double r) const{
        return exp(-r)/r;
    }
};

double circ(vector<double> x1, vector<double> x2){  //Umfang=Circumfrance
    return sqrt(pow(x1[0],2)+pow(x1[1],2))+sqrt(pow(x2[0],2)+pow(x2[1],2));
}

vector<double> sumvec(vector<double> x1, vector<double> x2, vector<double> x3=veczero){
    vector<double> x;

    for(int i=0; i<x3.size();++i) x.push_back(x1.at(i)+x2.at(i)+x3.at(i));
    return x;
}

vector<double> subvec(vector<double> x1, vector<double> x2){
    vector<double> x;
    for(int i=0; i<x1.size();++i) x.push_back(x1.at(i)-x2.at(i));
    return x;
}

vector<double> scalevec(vector<double> x1,double s){
    
    vector<double> x;
    for(int i=0; i<x1.size();++i) x.push_back(x1.at(i)*s);
    return x;
}

double euclid(vector<double> x){
    
    double sum=0;
    for(int i=0;i<x.size();i++){
        sum+=pow(x[i],2);
    }
    
    
    return sqrt(sum);
}




class LatticeSum{

    

     vector<vector<double>> minimizeCell(vector<double> x1, vector<double> x2) const{
        bool minimal=0;
        double u=circ(x1,x2);
        while(!minimal){
            bool xchange=0;
            double cells[4]={circ(sumvec(x1,x2,veczero),x2),circ(subvec(x1,x2),x2),circ(x1,sumvec(x2,x1,veczero)),circ(x1,subvec(x2,x1))};

            int smallest=-1;
            
            for(int i=0;i<4;i++){
                if(u>cells[i]){         
                    smallest=i;
                    u=cells[i];                        
                }
            }
             
            switch (smallest){
                case 0 : x1=sumvec(x1,x2);
                break;
                case 1 : x1=subvec(x1,x2);
                break;
                case 2 : x2=sumvec(x2,x1);
                break;
                case 3 : x2=subvec(x2,x1);
                break;
                case -1 : minimal=1;
                break;
                default : cout<<"Switch not found"<<endl;
            }
            u=circ(x1,x2);
        }
        return {x1,x2};
    }
    


    public:
    double operator() (double x, double phi, Potential v,bool min, vector<double> cx, vector<double> cy,vector<double> h) const{
        double sum=0;
        double a=1./sqrt(x*sin(phi)*density);
        vector<double> x1={a,0,0};
        vector<double> x2={a*x*cos(phi),a*x*sin(phi),0};
        vector<vector<double>> x3;

        for(int i=0;i<h.size();i++){
            vector<double> vec={cx[i],cy[i],h[i]};
            x3.push_back(vec);
        }

        /*
        if(min){
            vector<vector<double>> cell=minimizeCell(x1,x2);   
            x1=cell[0];
            x2=cell[1];
            x1.push_back(0);
            x2.push_back(0);
        }
        */
        
        int lmax=rcut/a/x;///sin(phi);
        lmax++;
        int kmax=rcut/a;
        kmax++;

        for(int k=-kmax;k<=kmax;++k){
                for(int l=-lmax;l<=lmax;++l){
                    if(!(l==0&&k==0)) sum+=v(euclid(sumvec(scalevec(x1,k),scalevec(x2,l))));
            }
        }
        sum*=layers;

        for(int i=0;i<=h.size();i++){
            //calculate layer
            vector<double> vertVec=veczero;
            for(int j=i-1;j>=0;j--){    //calculate layers under the current layer
                vertVec=sumvec(x3[j],vertVec);
                for(int k=-kmax;k<=kmax;++k){
                for(int l=-lmax;l<=lmax;++l){
                    sum+=v(euclid(sumvec(scalevec(x1,k),scalevec(x2,l),vertVec)));
                }}
            }
            vertVec=veczero;
            for(int j=i;j<x3.size();j++){    //calculate layers above the current layer
                vertVec=sumvec(x3[j],vertVec);
                for(int k=-kmax;k<=kmax;++k){
                for(int l=-lmax;l<=lmax;++l){
                    sum+=v(euclid(sumvec(scalevec(x1,k),scalevec(x2,l),vertVec)));
                }}
            }
        }
    
        return sum;
    }


     double operator() (double x, double phi, Potential v, vector<double> cx, vector<double> cy,vector<double> h) const{    //calculates just interaction energy between layers
        double sum=0;
        double a=1./sqrt(x*sin(phi)*density);
        vector<double> x1={a,0,0};
        vector<double> x2={a*x*cos(phi),a*x*sin(phi),0};
        vector<vector<double>> x3;

        for(int i=0;i<h.size();i++){
            vector<double> vec={cx[i],cy[i],h[i]};
            x3.push_back(vec);
        }

        /*
        if(min){
            vector<vector<double>> cell=minimizeCell(x1,x2);   
            x1=cell[0];
            x2=cell[1];
            x1.push_back(0);
            x2.push_back(0);
        }
        */
        
        int lmax=rcut/a/x;///sin(phi);
        lmax++;
        int kmax=rcut/a;
        kmax++;

        

        for(int i=0;i<=h.size();i++){
            //calculate layer
            vector<double> vertVec=veczero;
            for(int j=i-1;j>=0;j--){    //calculate layers under the current layer
                vertVec=sumvec(x3[j],vertVec);
                for(int k=-kmax;k<=kmax;++k){
                for(int l=-lmax;l<=lmax;++l){
                    sum+=v(euclid(sumvec(scalevec(x1,k),scalevec(x2,l),vertVec)));
                }}
            }
            vertVec=veczero;
            for(int j=i;j<x3.size();j++){    //calculate layers above the current layer
                vertVec=sumvec(x3[j],vertVec);
                for(int k=-kmax;k<=kmax;++k){
                for(int l=-lmax;l<=lmax;++l){
                    sum+=v(euclid(sumvec(scalevec(x1,k),scalevec(x2,l),vertVec)));
                }}
            }
        }
    
        return sum/layers;
    }
};

class Energy{

    LatticeSum latticeSum;
    Potential potential;
    public:
    double operator() (double x,double phi,vector<double> cx,vector<double> cy,vector<double> h=veczero){
        double sum=latticeSum(x,phi,potential,false,cx,cy,h);
        return sum/layers;
    }

};

class Fitness{

    
    LatticeSum latticeSum;
    Potential potential;
    public:
    double tric_sum;
    


        Fitness(){
            tric_sum=latticeSum(tric.x,tric.phi,potential,false,tric.cx,tric.cy,tric.hTric);
        }
        
        void calcTric(){
            tric_sum=latticeSum(tric.x,tric.phi,potential,false,tric.cx,tric.cy,tric.hTric);
        }
        
        double operator() (double x,double phi,bool min,vector<double> cx,vector<double> cy,vector<double> h=veczero) const{
            //if(phi<PI/6) return 0;  //deletes non sensical solutions
            double sum=latticeSum(x,phi,potential,min,cx,cy,h);
            return exp(1-pow(sum/tric_sum,1.+gen_i*0.05));
        }
} fit;





class Individual{
    bitset<bitlen> x;
    bitset<bitlen> phi;
    vector<bitset<bitlen>> cx;
    vector<bitset<bitlen>> cy;
    vector<bitset<bitlen>> h;

    double fitness=0;

    public:

    bitset<bitlen> getVecX(){
        return x;
    }
    bitset<bitlen> getVecPhi(){
        return phi;
    }

    vector<bitset<bitlen>> getBitCx(){
        return cx;
    }

    vector<bitset<bitlen>> getBitCy(){
        return cy;
    }

    vector<bitset<bitlen>> getBitH(){
        return h;
    }


    void setVecX(double dx){
        dx*=pow(2,x.size());
        long ix=--dx;

        bitset<bitlen> vecX;

        for(int i=0;i<x.size();i++){
            if(ix){
                vecX.set(i,ix&1);
                ix>>=1;
            }
            else vecX.set(i,0); 
        }
        
        //reverse(vecX.begin(),vecX.end());
        x=vecX;
    }

    void setVecPhi(double dx){
        dx/=(PI/2.);
        //if(dx>1) dx=1;
        dx*=pow(2,phi.size());
        long ix=--dx;
        
        bitset<bitlen> vecX;

        for(int i=0;i<phi.size();i++){
            if(ix){
                vecX.set(i,ix&1);
                ix>>=1;
            }
            else vecX.set(i,0); 
        }
        
        //reverse(vecX.begin(),vecX.end());
        phi=vecX;
    }

    void setBitCx(vector<double> dxVec){
        for(int i=0;i<dxVec.size();i++){
            double dx=dxVec[i];
            dx*=pow(2,cx[i].size());
            long ix=dx;

            bitset<bitlen> bitCx;

            for(int i=0;i<cx[i].size();++i){
                if(ix){
                    bitCx.set(i,ix&1);
                    ix>>=1;
                }
                else bitCx.set(i,0);
            }
            cx[i]=bitCx;
        }
    }

    void setBitCy(vector<double> dxVec){
        for(int i=0;i<dxVec.size();++i){
            double dx=dxVec[i];
            dx*=pow(2,cy[i].size());
            long ix=dx;

            bitset<bitlen> bitCy;

            for(int i=0;i<cy[i].size();++i){
                if(ix){
                    bitCy.set(i,ix&1);
                    ix>>=1;
                }
                else bitCy.set(i,0);
            }

            cy[i]=bitCy;
        }
    }

    void setBitH(vector<double> dhVec){
        for(int i=0;i<dhVec.size();i++){
            dhVec[i]*=pow(2,h[i].size());
            long iH=dhVec[i];

            bitset<bitlen> bitH;

            for(int j=0;j<h[i].size();j++){
                if(iH){
                    bitH.set(j,iH&1);
                    iH>>=1;
                }
                else bitH.set(i,0);
            }
            h[i]=bitH;
        }
    }



    double getX(){
        long accum=x.to_ulong();
        accum++;


        double ret=accum;
        ret/=pow(2,x.size());
        
        return ret;
    }
    double getPhi(){
        long accum=phi.to_ulong();
        accum++;

        double ret=accum;
        ret/=pow(2,phi.size());
        ret*=(PI/2);

        return ret;
    }

    double getA(){
        return 1./sqrt(getX()*sin(getPhi())*density);
    }


    //not sure if right
    vector<double> getCx(){
        vector<double> ret;
        for(int i=0;i<cx.size();i++){
            long accum=cx[i].to_ulong();

            
            double dret=accum;
            dret/=pow(2,cx[i].size());

            dret*=getA();

            double sub=0;
            if(i>0) for(int j=0;j<ret.size();j++) sub+=ret[j];
            dret-=sub;
            ret.push_back(dret);
        }
        return ret;
    }
    
    vector<double> getCy(){
        vector<double> ret;

        for(int i=0;i<cy.size();i++){
            long accum=cy[i].to_ulong();
        

            double dret=accum;
            dret/=pow(2,cy[i].size());
            dret*=getA()*getX()*sin(getPhi());
            
            double sub=0;
            if(i>0) for(int j=0;j<ret.size();j++) sub+=ret[j];
            dret-=sub;

            ret.push_back(dret);
        }
        return ret;
    }


    vector<double> getH(){
        vector<double> ret;


        for(int i=0;i<h.size();i++){
            long accum=h[i].to_ulong();

            double element=accum;
            element/=pow(2,h[i].size());
            ret.push_back(element);
        }

        double sum=0;
        for(int i=0;i<h.size();i++){
            sum+=ret[i];
        }
        for(int i=0;i<h.size();i++){
            ret[i]/=sum/d;
        }



        return ret;
    }


    double getFitness(){  
        return fitness;
    }

    void calcFitness(){
        fitness=fit(getX(),getPhi(),false,getCx(),getCy(),getH());
    }

    void printStats(){
        cout<<"Fitness="<<fitness<<endl<<"x="<<getX()<<endl<<"phi="<<getPhi()/PI*180<<endl<<"a="<<1./sqrt(sin(getPhi())*getX()*density)<<endl;
        cout<<"H values:"<<endl;
        vector<double> hVec=getH();
        vector<double> cxVec=getCx();
        vector<double> cyVec=getCy();
        for(int i=0;i<hVec.size();i++){
            cout<<"h="<<hVec[i]<<"; cx="<<cxVec[i]<<"; cy="<<cyVec[i]<<endl<<endl;
        }


        /*cout<<"X Genome=";
        for(int i=0;i<x.size();i++){
            cout<<x.at(i);
        }
        cout<<endl<<"Phi Genome=";
        for(int i=0;i<phi.size();i++){
            cout<<phi.at(i);
        }*/
        cout<<endl;
    }

    //Inital Creator
    Individual(int x_length,int phi_length){
        random_device rd;
        mt19937 gen(rd());
        uniform_int_distribution<> dis(0,1);

        fitness=1;
        for(int i=0;i<x.size();i++) x[i]=dis(gen);
        for(int i=0;i<phi.size();i++) phi[i]=dis(gen);
        
        for(int i=0;i<layers-1;i++){
            bitset<bitlen> bitH;
            bitset<bitlen> bitcx;
            bitset<bitlen> bitcy;

            for(int j=0;j<bitH.size();j++) bitH[j]=dis(gen);
            for(int i=0;i<bitcx.size();i++) bitcx[i]=dis(gen);
            for(int i=0;i<bitcy.size();i++) bitcy[i]=dis(gen);

            h.push_back(bitH);
            cx.push_back(bitcx);
            cy.push_back(bitcy);
        }
        correctCell();
    }
    //child Creator
    Individual(bitset<bitlen> xVec, bitset<bitlen> phiVec, vector<bitset<bitlen>> cx, vector<bitset<bitlen>> cy, vector<bitset<bitlen>> h){ //has to be changed for more or less than 2 parents!!!!!!
        x=xVec;
        phi=phiVec;
        this->cx=cx;
        this->cy=cy;
        this->h=h;
        mutate();

        correctCell();

    }
    Individual(){

        cout<<"I am called"<<endl;
        random_device rd;
        mt19937 gen(rd());
        uniform_int_distribution<> dis(0,1);

        fitness=1;
        for(int i=0;i<x.size();i++) x[i]=dis(gen);
        for(int i=0;i<phi.size();i++) phi[i]=dis(gen);
        for(int i=0;i<cx.size();i++) cx[i]=dis(gen);
        for(int i=0;i<cy.size();i++) cy[i]=dis(gen);
        for(int i=0;i<layers-1;i++){
            bitset<bitlen> bit;
            for(int j=0;j<bit.size()-1;j++) bit[j]=dis(gen);
            h.push_back(bit);
        }
        correctCell();
    }

    private:
    void mutate(){
        random_device rd;
        mt19937 gen(rd());
        uniform_real_distribution<> dis(0,1);
        for(int i=0;i<x.size();++i){
            if(dis(gen)<=mutation_rate){
                x.flip(i);
            }
        }
        for(int i=0;i<phi.size();++i){
            if(dis(gen)<=mutation_rate){
                phi.flip(i);
            }
        }
        for(int j=0;j<cx.size();++j){
            for(int i=0;i<cx[j].size();++i){
                if(dis(gen)<=mutation_rate){
                    cx[j].flip(i);
                }
            }
        }
        
        for(int j=0;j<cy.size();++j){
            for(int i=0;i<cy[j].size();++i){
                if(dis(gen)<=mutation_rate){
                    cy[j].flip(i);
                }
            }
        }

        for(int i=0;i<h.size();++i){
            for(int j=0;j<h[i].size();++j){
                if(dis(gen)<=mutation_rate){
                    h[i].flip(j);
                }
            }
        }


    }
    public: 
    bool corrected=false;
    double lol=90;
    double tol=2;

public:
    void correctCell(){
        double x=getX();
        double phi=getPhi();
        

        double a=1./sqrt(x*sin(phi)*density);
        vector<double> x1={a,0,0};
        vector<double> x2={a*x*cos(phi),a*x*sin(phi),0};
        
        
        vector<vector<double>> vecX=minimizeCell(x1,x2);
        x1=vecX[0];
        x2=vecX[1];

        if(euclid(x1)<euclid(x2)){
            //cout<<phi/PI*180<<endl;
            vector<double> temp=x1;
            x1=x2;
            x2=temp;
        }
        
        
        a=euclid(x1);
        //cout<<a<<endl;
        x=euclid(x2)/a;
        
        phi=acos((x1[0]*x2[0]+x1[1]*x2[1])/(euclid(x1)*euclid(x2)));
        if(phi>PI/2){
            x2=sumvec(x2,x1);
            x=euclid(x2)/a;
            phi=acos((x1[0]*x2[0]+x1[1]*x2[1])/(euclid(x1)*euclid(x2)));
        } 

        vector<double> cx=getCx();
        vector<double> cy=getCy();

        for(int i=0;i<cx.size();i++){
            /*if(cx[i]>a) cx[i]-=a;
            if(cy[i]>a) cy[i]-=a;*/
            if(abs(cx[i])<abs(cy[i])){
                bitset<bitlen> temp=(this->cx)[i];
                (this->cx)[i]=(this->cy)[i];
                (this->cy)[i]=temp;
            }
            
        }
        
        vector<double> h=getH();

        if(h[0]<h[h.size()-1]){
            reverse((this->h).begin(),(this->h).end());
            reverse((this->cx).begin(),(this->cx).end());
            reverse((this->cy).begin(),(this->cy).end());
        }


        setVecX(x);
        setVecPhi(phi);

        corrected=true;
        lol=phi;
    }




    vector<vector<double>> minimizeCell(vector<double> x1, vector<double> x2) const{
        bool minimal=0;
        double u=circ(x1,x2);
        while(!minimal){
            /*cout<<x1[0]<<";"<<x1[1]<<endl;
            cout<<x2[0]<<";"<<x2[1]<<endl<<endl;*/
            double cells[4]={circ(sumvec(x1,x2),x2),circ(subvec(x1,x2),x2),circ(x1,sumvec(x2,x1)),circ(x1,subvec(x2,x1))};

            int smallest=-1;
            
            for(int i=0;i<4;i++){
                if(u>cells[i]){         
                    smallest=i;
                    u=cells[i];                        
                }
            }
             
            switch (smallest){
                case 0 : x1=sumvec(x1,x2);
                break;
                case 1 : x1=subvec(x1,x2);
                break;
                case 2 : x2=sumvec(x2,x1);
                break;
                case 3 : x2=subvec(x2,x1);
                break;
                case -1 : minimal=1;
                break;
                default : cout<<"Switch not found"<<endl;
            }
            u=circ(x1,x2);
        }
        return {x1,x2};
    }


};






bool fitComp(Individual a, Individual b){ return a.getFitness()<b.getFitness();}




class Generation{
    vector<Individual> indiv;
    double genFitness=0;

    mutex child_lock;

    void createChildren(vector<Individual> *children,int size){
        for(int i=0;i<size;i+=2){      // double because 2 children are created
            vector<Individual> child=genChild();
            //cout<<"First child generated"<<endl;
            child_lock.lock();
            children->insert(children->end(),child.begin(),child.end());
            child_lock.unlock();
        }
    }

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


    Individual getBest(){
        Individual best=indiv.at(0);
        for(int i=1;i<indiv.size();i++){
            if(indiv.at(i).getFitness()>best.getFitness()) best=indiv.at(i);
        }
        return best;
    }

    void printBest(){
        Individual best=getBest();
        best.printStats();
    }


    void nextGen(int size){
       calcGenFitness();
        //printBest();
       vector<Individual> children;
        //cout<<"children will now be born"<<endl;
        

        vector<thread> t;
        if(threadcount==1){
            createChildren(&children,(size-survivors));
        }
        else{
        for(int i=0;i<threadcount;++i){
            t.push_back(thread(&Generation::createChildren,this,&children,(size-survivors)/threadcount));       //size has to be devisable by threadcoutn
        }
        for(int i=0;i<threadcount;++i){
            t.at(i).join();
        }
        }

        vector<Individual> surv=getSurvivors(survivors);
        children.insert(children.end(),surv.begin(),surv.end());
        
        //cout<<"Children become Parents"<<endl;
        indiv=children;
        children.clear();
        //printBest();
    }
    
    vector<Individual> getSurvivors(int count){
        vector<Individual> survivors=indiv;
        std::sort(indiv.begin(),indiv.end(),fitComp);
        survivors.resize(count);
        return survivors;
    }


    
    /*
    void printIndiv(){
        calcGenFitness();
        auto print=[](Individual ind){cout<<"Individual:"<<endl<<"X="<<ind.getX()<<endl<<"Phi="<<ind.getPhi()/PI*180<<endl<<"Cx="<<ind.getCx()<<endl<<"Cy="<<ind.getCy()<<endl<<"Fit="<<ind.getFitness()<<endl<<endl;};
        for_each(indiv.begin(),indiv.end(),print);
    }
*/
    Individual getLastBest(){
        calcGenFitness();
        printBest();
        return getBest();
    }

    

    private:

    Individual getParent(){
        
        random_device rd;
        mt19937 gen(rd());
        uniform_real_distribution<> dis(0,genFitness);     
        double rand=dis(gen);
        int i=0;
        
        
        
        while(indiv.at(i).getFitness()<rand){
            //cout<<i<<" "<<indiv.at(i).getFitness()<<endl;
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
        uniform_int_distribution<> disCx(0,bitlen);

        int cutX=disX(gen);
        int cutPhi=disPhi(gen);

        vector<int> cutCx;
        for(int i=0;i<parents.at(0).getBitCx().size();++i){
            cutCx.push_back(disCx(gen));
        }

        vector<int> cutCy;
        for(int i=0;i<parents.at(0).getBitCy().size();++i){
            cutCy.push_back(disCx(gen));
        }
        
        vector<int> cutH;
        for(int i=0;i<parents.at(0).getBitH().size();++i){
            cutH.push_back(disCx(gen));
        }

        
        bitset<bitlen> x1;
        bitset<bitlen> phi1;
        bitset<bitlen> x2;
        bitset<bitlen> phi2;
        
        vector<bitset<bitlen>> cx1;
        vector<bitset<bitlen>> cx2;
        vector<bitset<bitlen>> cy1;
        vector<bitset<bitlen>> cy2;

        vector<bitset<bitlen>> h1;
        vector<bitset<bitlen>> h2;
        
        /*using namespace std::chrono;
        high_resolution_clock::time_point t_start_parallel = high_resolution_clock::now();
        *///cout<<"Genes will be created!"<<endl;
        for(int i=0;i<cutX;++i){
            x1[i]=parents[0].getVecX()[i];
            x2[i]=parents[1].getVecX()[i];
        }
        //cout<<"Secon Gene now edited"<<endl;
        for(int i=cutX;i<parents.at(0).getVecX().size();++i){
            x1[i]=parents[1].getVecX()[i];
            x2[i]=parents[0].getVecX()[i];
        }

        for(int i=0;i<cutPhi;++i){
            phi1[i]=parents[0].getVecPhi()[i];
            phi2[i]=parents[1].getVecPhi()[i];
        }
        for(int i=cutPhi;i<parents.at(0).getVecPhi().size();++i){
            phi1[i]=parents[1].getVecPhi()[i];
            phi2[i]=parents[0].getVecPhi()[i];
        }


        for(int j=0;j<cutCx.size();++j){
            bitset<bitlen> bit1;
            bitset<bitlen> bit2;


            for(int i=0;i<cutCx[j];++i){
                bit1[i]=parents[0].getBitCx()[j][i];
                bit2[i]=parents[1].getBitCx()[j][i];
            }
            for(int i=cutCx[j];i<bitlen;++i){
                bit1[i]=parents[0].getBitCx()[j][i];
                bit2[i]=parents[1].getBitCx()[j][i];
            }
            cx1.push_back(bit1);
            cx2.push_back(bit2);
        }
        for(int j=0;j<cutCy.size();++j){
            bitset<bitlen> bit1;
            bitset<bitlen> bit2;
            for(int i=0;i<cutCy[j];++i){
                bit1[i]=parents[0].getBitCy()[j][i];
                bit2[i]=parents[1].getBitCy()[j][i];
            }
            for(int i=cutCy[j];i<bitlen;++i){
                bit1[i]=parents[0].getBitCy()[j][i];
                bit2[i]=parents[1].getBitCy()[j][i];
            }
            cy1.push_back(bit1);
            cy2.push_back(bit2);
        }




        for(int i=0;i<cutH.size();++i){
            bitset<bitlen> bit1;
            bitset<bitlen> bit2;


            for(int j=0;j<cutH[i];++j){
                bit1[j]=parents[0].getBitH()[i][j];
                bit2[j]=parents[1].getBitH()[i][j];
            }
            for(int j=cutH[i];j<bitlen;++j){
                bit2[j]=parents[0].getBitH()[i][j];
                bit1[j]=parents[1].getBitH()[i][j];
            }
            h1.push_back(bit1);
            h2.push_back(bit2);
        }

        
        /*high_resolution_clock::time_point t_end_parallel = high_resolution_clock::now();
        duration<double> time_parallel = t_end_parallel - t_start_parallel;
        cout << "Execution time: " << time_parallel.count()<<endl;
        */

        //cout<<"The two children will be created"<<endl;
        Individual child1(x1,phi1,cx1,cy1,h1);
        Individual child2(x2,phi2,cx2,cy2,h2);

        vector<Individual> children;
        children.push_back(child1);
        children.push_back(child2);
       
        return children;
    }

    mutex indexMutex;

    void threadCalcFitness(vector<Individual> *indiv, int *index){
        while(true){
            indexMutex.lock();
            int i=(*index)--;
            indexMutex.unlock();

            if(i<0) break;

            indiv->at(i).calcFitness();

        }
    }

    void calcGenFitness(){
        int nextIndex=indiv.size()-1;

        vector<thread> t;

        for(int i=0;i<threadcount;i++){
            t.push_back(thread(&Generation::threadCalcFitness,this,&indiv,&nextIndex));
        }
        for(int i=0;i<threadcount;i++){
            t.at(i).join();
        }


        genFitness=accumulate(indiv.begin(),indiv.end(),0.,
            [](double x, Individual y){
            //y.calcFitness(); 
            //cout<<y.getFitness()<<endl;
            return x+y.getFitness();
            });
    }

    

};


class Climber{
    double x;
    double phi;
    Fitness fit;

    public:

    vector<double> hillclimb(double x,double phi, Fitness fit, vector<double> cx, vector<double> cy, vector<double> h){
        double stepX=0.001;
        double stepPhi=0.0001;
        
        bool top=0;
        double bestfit=0;

        while(!top){

            
            



            while(fit(x+stepX,phi,true,cx,cy,h)>fit(x,phi,true,cx,cy,h)&&x+stepX<=1) x+=stepX;
            while(fit(x-stepX,phi,true,cx,cy,h)>fit(x,phi,true,cx,cy,h)&&x-stepX>=0) x-=stepX;
            
            /*
            while(fit(x,phi,true,cx+stepX,cy,h)>fit(x,phi,true,cx,cy,h)) cx+=stepX;
            while(fit(x,phi,true,cx-stepX,cy,h)>fit(x,phi,true,cx,cy,h)) cx-=stepX;
            while(fit(x,phi,true,cx,cy+stepX,h)>fit(x,phi,true,cx,cy,h)) cy+=stepX;
            while(fit(x,phi,true,cx,cy-stepX,h)>fit(x,phi,true,cx,cy,h)) cy-=stepX;
            */


            while(fit(x,phi+stepPhi,true,cx,cy,h)>fit(x,phi,true,cx,cy,h)&&phi+stepPhi<=PI/2) phi+=stepPhi;
            while(fit(x,phi-stepPhi,true,cx,cy,h)>fit(x,phi,true,cx,cy,h)&&phi-stepPhi>0) phi-=stepPhi;

            double x1=1./sqrt(x*sin(phi)*density);
            double x2=x1*x*sin(phi);


            for(int i=0;i<cx.size();i++){
                vector<double> cxPlus=cx;
                cxPlus[i]+=stepX;

                double adder=0; //for cx addition of other layers
                for(double x:cxPlus) adder+=x;
                while(fit(x,phi,true,cxPlus,cy,h)>fit(x,phi,true,cx,cy,h)&&(adder<=x1)&&cxPlus[0]<x1){
                    cxPlus[i]+=stepX;
                    cx[i]+=stepX;
                    
                    adder=0;
                    for(double x:cxPlus) adder+=x;
                }

                vector<double> cxMinus=cx;
                
                adder=0;
                for(double x:cxMinus) adder+=x;
                cxMinus[i]-=stepX;
                while(fit(x,phi,true,cxMinus,cy,h)>fit(x,phi,true,cx,cy,h)&&adder>=0&&cxMinus[0]>0){
                    cxMinus[i]-=stepX;
                    cx[i]-=stepX;
                    adder=0;
                    for(double x:cxMinus) adder+=x;
                }

                
                vector<double> cyPlus=cy;
                adder=0;
                for(double x:cyPlus) adder+=x;
                cyPlus[i]+=stepX;
                while(fit(x,phi,true,cx,cyPlus,h)>fit(x,phi,true,cx,cy,h)&&adder<=x2&&cxPlus[0]<x2){
                    cyPlus[i]+=stepX;
                    cy[i]+=stepX;
                    adder=0;
                    for(double x:cyPlus) adder+=x;
                }

                
                vector<double> cyMinus=cy;
                adder=0;
                for(double x:cyMinus) adder+=x;
                cyMinus[i]-=stepX;
                while(fit(x,phi,true,cx,cyMinus,h)>fit(x,phi,true,cx,cy,h)&&adder>=0&&cyMinus[0]>0){
                    cyMinus[i]-=stepX;
                    cy[i]-=stepX;
                    adder=0;
                    for(double x:cyMinus) adder+=x;
                }


            }


            for(int i=0;i<h.size();i++){
                vector<double> hplus=h;
                hplus[i]+=stepX;

                double sum=0;
                for(int i=0;i<h.size();i++){
                    sum+=hplus[i];
                }
                for(int i=0;i<h.size();i++){
                    hplus[i]/=sum/d;
                }
                sum=0;

                
                while(fit(x,phi,true,cx,cy,hplus)>fit(x,phi,true,cx,cy,h)){
                    h=hplus;
                    hplus[i]+=stepX;

                    for(int i=0;i<h.size();i++){
                        sum+=hplus[i];
                    }
                    for(int i=0;i<h.size();i++){
                        hplus[i]/=sum/d;
                    }
                    sum=0;

                } 

                vector<double> hminus=h;
                hminus[i]-=stepX;
                sum=0;
                for(int i=0;i<h.size();i++){
                    sum+=hminus[i];
                }
                for(int i=0;i<h.size();i++){
                    hminus[i]/=sum/d;
                }
                sum=0;


                while(fit(x,phi,true,cx,cy,hminus)>fit(x,phi,true,cx,cy,h)&&hminus[i]>=0){
                    h=hminus;
                    hminus[i]-=stepX;

                    for(int i=0;i<h.size();i++){
                        sum+=hminus[i];
                    }
                    for(int i=0;i<h.size();i++){
                        hminus[i]/=sum/d;
                    }
                    sum=0;

                } 


            }


            if(bestfit-fit(x,phi,true,cx,cy,h)<0.0001) top=1;
            bestfit=fit(x,phi,true,cx,cy,h);
        }
        vector<double> ret={x,phi};
        for(int i=0;i<layers-1;i++){
            ret.push_back(h[i]);
            ret.push_back(cx[i]);
            ret.push_back(cy[i]);
        }
        return ret;
    }
};




namespace plt=matplotlibcpp;

void plotCell(vector<double> top,string plotname="crystall.png",double energy=0){
    double x=top[0];
    double phi=top[1];
    top.erase(top.begin(),top.begin()+2);


    double a=1./sqrt(x*sin(phi)*density);
    vector<double> x1={a,0,0};
    vector<double> x2={a*x*cos(phi),a*x*sin(phi),0};
    
    vector<vector<double>> c;

    for(int i=0;i<top.size();i+=3){
        c.push_back({top[i+1],top[i+2],top[i]});
    }

    vector<double> p1,p2,p3;
    vector<double> u1,u2,u3;
    vector<double> up1,up2,up3;
    vector<double> upp1,upp2,upp3;
    vector<double> p5_1,p5_2,p5_3;
    vector<double> p6_1,p6_2,p6_3;

    
    int lmax=3;
    int kmax=3;


    plt::clf();

    for(int i=-layers/2;i<=(layers+1)/2;i++){
    for(int k=-kmax;k<=kmax;++k){
        for(int l=-lmax;l<=lmax;++l){

            int summer=layers-1;

            switch (layers){
                case 6:{
                    vector<double> sumC=c[0];
                    for(int i=1;i<summer;i++) sumC=sumvec(sumC,c[i]);

                    vector<double> uppper_vec=sumvec(scalevec(x1,k),scalevec(x2,l),sumC);
                    p6_1.push_back(uppper_vec[0]);
                    p6_2.push_back(uppper_vec[1]);
                    p6_3.push_back(uppper_vec[2]);

                    summer--;
                }

                case 5:{
                    vector<double> sumC=c[0];
                    for(int i=1;i<summer;i++) sumC=sumvec(sumC,c[i]);

                    vector<double> uppper_vec=sumvec(scalevec(x1,k),scalevec(x2,l),sumC);
                    p5_1.push_back(uppper_vec[0]);
                    p5_2.push_back(uppper_vec[1]);
                    p5_3.push_back(uppper_vec[2]);

                    summer--;
                }

                case 4:{
                    vector<double> sumC=c[0];
                    for(int i=1;i<summer;i++) sumC=sumvec(sumC,c[i]);

                    vector<double> uppper_vec=sumvec(scalevec(x1,k),scalevec(x2,l),sumC);
                    upp1.push_back(uppper_vec[0]);
                    upp2.push_back(uppper_vec[1]);
                    upp3.push_back(uppper_vec[2]);

                    summer--;
                }
                case 3: {
                    vector<double> sumC=c[0];
                    for(int i=1;i<summer;i++) sumC=sumvec(sumC,c[i]);

                    vector<double> upper_vec=sumvec(scalevec(x1,k),scalevec(x2,l),sumC);
                    up1.push_back(upper_vec[0]);
                    up2.push_back(upper_vec[1]);
                    up3.push_back(upper_vec[2]);
                    summer--;
                }
                case 2: {
                    vector<double> sumC=c[0];


                    vector<double> under_vec=sumvec(scalevec(x1,k),scalevec(x2,l),sumC);
                    u1.push_back(under_vec[0]);
                    u2.push_back(under_vec[1]);
                    u3.push_back(under_vec[2]);


                }
                case 1: {


                    vector<double> vec=sumvec(scalevec(x1,k),scalevec(x2,l));
                    p1.push_back(vec[0]);
                    p2.push_back(vec[1]);
                    p3.push_back(0);
                }
                default: break;
            }
        }
    }
    }
   
    plt::plot(p1,p2,"ro",u1,u2,"bs",up1,up2,"gv",upp1,upp2,"yh",p5_1,p5_2,"rs",p6_1,p6_2,"bo");
    //string title="x="+to_string(x)+" , phi="+to_string(phi/PI*180.)+" , fitness="+to_string(fitness);
    string title="d="+to_string(d)+", l="+to_string(layers)+", vD="+to_string(volDensity)+", lD="+to_string(density)+", E="+to_string(energy);
    plt::title(title);
    //cout<<bo;
    //bo=plt::plot(p1,p2);
    //cout<<bo;
    plt::axis("equal");
    //plt::show();

    string plotname1=plotname+"_1.png";
    plt::save(plotname1);

    plt::clf();

    plt::plot(p1,p3,"ro",u1,u3,"bs",up1,up3,"gv",upp1,upp3,"yh",p5_1,p5_3,"rs",p6_1,p6_3,"bo");
    plt::title(title);

    string plotname2=plotname+"_2.png";
    plt::save(plotname2);
}



int main(){

    const int ind_size=4000;
    const int generations=300;

    for(volDensity=0.1;volDensity<=0.1;volDensity+=0.02){
        for(d=3.6;d<=3.8;d+=0.01){
            for(layers=2;layers<=6;layers+=1){
                cout<<"Now calculating dens="<<volDensity<<" d="<<d<<" layers="<<layers<<endl;

                density=d*volDensity/layers;
                tric.calculateTric();
                fit.calcTric();

                Generation gen(ind_size,8,6);

                for(int i=0;i<generations;i++){
                    gen_i=i+1;  //for better GA performance, see fitness function
                    if(!(i%10)) cout<<"Generation "<<i<<" of "<<generations<<endl<<endl;

                    gen.nextGen(ind_size);
                }
                
                cout<<"GA calculated stats:"<<endl;
                Individual best=gen.getLastBest();

                Climber climb;
                vector<double> top= climb.hillclimb(best.getX(),best.getPhi(),fit,best.getCx(),best.getCy(),best.getH());
                cout<<"After he climbed the hill:"<<endl<<"X="<<top[0]<<endl<<"Phi="<<top[1]*180/PI<<endl<<"a="<<(1./sqrt(density*top[0]*sin(top[1])))<<endl;
                cout<<"H Values(h,cx,cy):"<<endl;


                vector<double> hFit;
                vector<double> cxFit;
                vector<double> cyFit;
                for(int i=2;i<top.size();i+=3){
                    cout<<(i-2)/3<<":"<<endl;
                    cout<<top[i]<<endl;
                    cout<<top[i+1]<<endl;
                    cout<<top[i+2]<<endl<<endl;
                    hFit.push_back(top[i]);
                    cxFit.push_back(top[i+1]);
                    cyFit.push_back(top[i+2]);
                }
                cout<<"Fitness="<<fit(top[0],top[1],false,cxFit,cyFit,hFit)<<endl;
                Energy energy;
                double energy_value=energy(top[0],top[1],cxFit,cyFit,hFit);
                cout<<"Energy="<<energy_value<<endl;
                LatticeSum latticeSum;
                Potential v;
                cout<<"Energy of Interaction between Layers="<<latticeSum(top[0],top[1],v,cxFit,cyFit,hFit)<<endl<<endl;


                string plotname="dens="+to_string(volDensity)+"_d="+to_string(d)+"_l="+to_string(layers);
                plotCell(top,plotname,energy_value);
                
            }
        }
    }
/*
    Fitness fit;

    Generation gen(ind_size,8,6);

  
   
    for(int j=1;j<=1000;j++){
        /*d=0.1*j*lambda;
        tric.setH();
        using namespace std::chrono;
        high_resolution_clock::time_point t_start_parallel = high_resolution_clock::now();

        int imax=200;

        for(int i=0;i<imax;i++){
            gen_i=i+1;
            //cout<<gen_i<<" ";
            //d=j*0.1;
                //high_resolution_clock::time_point t_start_parallel = high_resolution_clock::now();
            cout<<endl<<"Generation "<<i<<" of "<<imax<<endl;
            gen.nextGen(ind_size);
                //high_resolution_clock::time_point t_end_parallel = high_resolution_clock::now();
                //duration<double> time_parallel = t_end_parallel - t_start_parallel;
                //cout << "Execution time: " << time_parallel.count()<<endl;
            
            
            Individual best=gen.getBest();
            //best.printStats();
            string plotname="Generation";
            plotname+=to_string(i);
            plotname+=".png";
            best.calcFitness();
                //plotCell(best.getX(),best.getPhi(),plotname,best.getFitness());
            
        }
        cout<<endl;
        //gen.printIndiv();
        cout<<endl<<endl<<"Winners of the Evolution Contest "<<j<<" :"<<endl;
        Individual best=gen.getLastBest();
        //best.correctCell();

        //best.printStats();
        //cout<<endl<<endl<<"ALL THE INDIVS"<<endl;
        //gen.printIndiv();
        Climber climb;
        vector<double> top= climb.hillclimb(best.getX(),best.getPhi(),fit,best.getCx(),best.getCy(),best.getH());
        cout<<endl<<endl<<"After he climbed the hill:"<<endl<<"X="<<top[0]<<endl<<"Phi="<<top[1]*180/PI<<endl<<"a="<<(1./sqrt(density*top[0]*sin(top[1])))<<endl;
        cout<<"H Values(h,cx,cy):"<<endl;

        vector<double> hFit;
        vector<double> cxFit;
        vector<double> cyFit;
        for(int i=2;i<top.size();i+=3){
            cout<<(i-2)/3<<":"<<endl;
            cout<<top[i]<<endl;
            cout<<top[i+1]<<endl;
            cout<<top[i+2]<<endl<<endl;
            hFit.push_back(top[i]);
            cxFit.push_back(top[i+1]);
            cyFit.push_back(top[i+2]);
        }
        cout<<"Fitness="<<fit(top[0],top[1],false,cxFit,cyFit,hFit)<<endl;
        Energy energy;
        double energy_value=energy(top[0],top[1],cxFit,cyFit,hFit);
        cout<<"Energy="<<energy_value<<endl;
        string plotname="crystall";
        plotname+=to_string(j);
        plotname+=".png";
        plotCell(top,plotname,energy_value);

        high_resolution_clock::time_point t_end_parallel = high_resolution_clock::now();
        duration<double> time_parallel = t_end_parallel - t_start_parallel;
        cout << "Execution time: " << time_parallel.count()<<endl;



        gen_i=1;
    }/*
    vector<bool> x={0,1,1,0,1,1,1,0};
    vector<bool> phi={0,0,0,0,0,0,0};
    Individual a(8,6);
    a.setVecPhi(PI/3);
    a.correctCell();
    cout<<a.getPhi()<<endl;*/
    
    
    cout<<"FINISHED! \nEnter anything to close"<<endl;
    bool close;
    cin>>close;
    return 0;
}