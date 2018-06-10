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
double d=lambda*0.05;
int n=100;
double gen_i=1;
int const bitlen=10;
const double density=1;

int const layers=2;

double rcut=lambda*13;
//int kmax=15;
//int lmax=15;

int parent_amount=2;        //changes in code have to be made!!!
double mutation_rate=0.1;

vector<double> veczero={0,0,0};

int threadcount=4;

int survivors=1;//this many best genes will deifnitly survive

struct Tric{
    public:
    double const phi=PI/2.;
    double const x=1;
    double const a=1./sqrt(x*sin(phi)*density);
    double const cx=a/2.;
    double const cy=a*x/2.;
    vector<double> hTric;

    Tric(){
        setH();
    }

    void setH(){
        hTric.erase(hTric.begin(),hTric.end());
        for(int i=0;i<layers-1;i++){
            hTric.push_back(d/(layers-1.));
        }
    }



} tric;

class Potential{
    public:
    double operator() (double r) const{
        return exp(-r/lambda)/r;
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
    double operator() (double x, double phi, Potential v,bool min=false, double cx=0, double cy=0,vector<double> h=veczero) const{
        double sum=0;
        double a=1./sqrt(x*sin(phi)*density);
        vector<double> x1={a,0,0};
        vector<double> x2={a*x*cos(phi),a*x*sin(phi),0};
        vector<vector<double>> x3;

        for(int i=h.size()-1;i>=0;i--){
            vector<double> vec={cx,cy,h[i]};
            x3.push_back(vec);
        }

        
        if(min){
            vector<vector<double>> cell=minimizeCell(x1,x2);   
            x1=cell[0];
            x2=cell[1];
            x1.push_back(0);
            x2.push_back(0);
        }
        
        
        int lmax=rcut/a/x;///sin(phi);
        lmax++;
        int kmax=rcut/a;
        kmax++;

        for(int i=0;i<h.size();i++){
            //calculate layer
            for(int k=-kmax;k<=kmax;++k){
            for(int l=-lmax;l<=lmax;++l){
                if(!(l==0&&k==0)) sum+=v(euclid(sumvec(scalevec(x1,k),scalevec(x2,l))));
            }
            }
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
};


class Fitness{

    
    LatticeSum latticeSum;
    Potential potential;
    public:
    double tric_sum;
    


        Fitness(){
            tric_sum=latticeSum(tric.x,tric.phi,potential,false,tric.cx,tric.cy,tric.hTric);
        }
        
        double operator() (double x,double phi,bool min=false,double cx=0,double cy=0,vector<double> h=veczero) const{
            //if(phi<PI/6) return 0;  //deletes non sensical solutions
            double sum=latticeSum(x,phi,potential,min,cx,cy,h);
            return exp(1-pow(sum/tric_sum,1.+gen_i*0.1));
        }
} fit;





class Individual{
    bitset<bitlen> x;
    bitset<bitlen> phi;
    bitset<bitlen> cx;
    bitset<bitlen> cy;
    vector<bitset<bitlen>> h;

    double fitness=0;

    public:

    bitset<bitlen> getVecX(){
        return x;
    }
    bitset<bitlen> getVecPhi(){
        return phi;
    }

    bitset<bitlen> getBitCx(){
        return cx;
    }

    bitset<bitlen> getBitCy(){
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
        dx*=pow(2,phi.size());
        dx/=(PI/2);
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

    void setBitCx(double dx){
        dx*=pow(2,cx.size());
        long ix=--dx;

        bitset<bitlen> bitCx;

        for(int i=0;i<cx.size();++i){
            if(ix){
                bitCx.set(i,ix&1);
                ix>>=1;
            }
            else bitCx.set(i,0);
        }

        cx=bitCx;
    }

    void setBitCy(double dx){
        dx*=pow(2,cy.size());
        long ix=--dx;

        bitset<bitlen> bitCy;

        for(int i=0;i<cy.size();++i){
            if(ix){
                bitCy.set(i,ix&1);
                ix>>=1;
            }
            else bitCy.set(i,0);
        }

        cy=bitCy;
    }

    void setBitH(vector<double> dhVec){
        for(int i=0;i<dhVec.size();i++){
            dhVec[i]*=pow(2,h[i].size());
            long iH=--dhVec[i];

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
        ret*=(PI/2);
        ret/=pow(2,phi.size());
        
        return ret;
    }


    //not sure if right
    double getCx(){
        long accum=cx.to_ulong();
        accum++;

        
        double ret=accum;
        ret/=pow(2,cx.size());
        return ret;
    }
    
    double getCy(){
        long accum=cy.to_ulong();
        accum++;
     

        double ret=accum;
        ret/=pow(2,cy.size());
        return ret;
    }


    vector<double> getH(){
        vector<double> ret;


        for(int i=h.size()-1;i>=0;i--){
            long accum=h[i].to_ulong();
            accum++;

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
        cout<<"Fitness="<<fitness<<endl<<"x="<<getX()<<endl<<"phi="<<getPhi()/PI*180<<endl<<"a="<<1./sqrt(sin(getPhi())*getX()*density)<<endl<<"cx="<<getCx()<<endl<<"cy="<<getCy()<<endl;
        cout<<"H values:"<<endl;
        vector<double> hVec=getH();
        for(int i=0;i<hVec.size();i++){
            cout<<hVec[i]<<endl;
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
        for(int i=0;i<cx.size();i++) cx[i]=dis(gen);
        for(int i=0;i<cy.size();i++) cy[i]=dis(gen);
        for(int i=0;i<layers-1;i++){
            bitset<bitlen> bit;
            for(int j=0;j<bit.size()-1;j++) bit[j]=dis(gen);
            h.push_back(bit);
        }
        correctCell();
    }
    //child Creator
    Individual(bitset<bitlen> xVec, bitset<bitlen> phiVec, bitset<bitlen> cx, bitset<bitlen> cy, vector<bitset<bitlen>> h){ //has to be changed for more or less than 2 parents!!!!!!
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
        for(int i=0;i<cx.size();++i){
            if(dis(gen)<=mutation_rate){
                cx.flip(i);
            }
        }
        for(int i=0;i<cy.size();++i){
            if(dis(gen)<=mutation_rate){
                cy.flip(i);
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
        
        phi=acos(x2[0]/(x*a));
        //cout<<phi<<endl;
        /*if(phi<PI/3){
            cout<<"Strange things are happening"<<endl;
            cout<<x1[0]<<";"<<x1[1]<<";"<<a<<endl;
            cout<<x2[0]<<";"<<x2[1]<<";"<<euclid(x2)<<endl<<endl;
        }*/

/*
        x1={a,0,0};
        x2={a*x*cos(phi),a*x*sin(phi),0};
        vecX=minimizeCell(x1,x2);
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
        phi=asin(1./(density*a*a*x));
*/





        double cx=getCx();
        double cy=getCy();
        if(cy>cx){
            bitset<bitlen> temp=this->cx;
            this->cx=this->cy;
            this->cy=temp;
        }

        setVecX(x);
        setVecPhi(phi);

        corrected=true;
        lol=phi;
        tol=density*a*a*x;
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

        for(int i=0;i<threadcount;++i){
            t.push_back(thread(&Generation::createChildren,this,&children,(size-survivors)/threadcount));       //size has to be devisable by threadcoutn
        }
        for(int i=0;i<threadcount;++i){
            t.at(i).join();
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


    

    void printIndiv(){
        calcGenFitness();
        auto print=[](Individual ind){cout<<"Individual:"<<endl<<"X="<<ind.getX()<<endl<<"Phi="<<ind.getPhi()/PI*180<<endl<<"Cx="<<ind.getCx()<<endl<<"Cy="<<ind.getCy()<<endl<<"Fit="<<ind.getFitness()<<endl<<endl;};
        for_each(indiv.begin(),indiv.end(),print);
    }

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
        int cutCx=disCx(gen);
        int cutCy=disCx(gen);
        
        vector<int> cutH;
        for(int i=0;i<parents.at(0).getBitH().size();++i){
            cutH.push_back(disCx(gen));
        }

        
        bitset<bitlen> x1(parents.at(0).getVecX().size());
        bitset<bitlen> phi1(parents.at(0).getVecPhi().size());
        bitset<bitlen> x2(parents.at(0).getVecX().size());
        bitset<bitlen> phi2(parents.at(0).getVecPhi().size());
        
        bitset<bitlen> cx1;
        bitset<bitlen> cx2;
        bitset<bitlen> cy1;
        bitset<bitlen> cy2;

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

        for(int i=0;i<cutCx;++i){
            cx1[i]=parents[0].getBitCx()[i];
            cx2[i]=parents[1].getBitCx()[i];
        }
        for(int i=cutCx;i<bitlen;++i){
            cx2[i]=parents[0].getBitCx()[i];
            cx1[i]=parents[1].getBitCx()[i];
        }
        
        for(int i=0;i<cutCy;++i){
            cy1[i]=parents[0].getBitCy()[i];
            cy2[i]=parents[1].getBitCy()[i];
        }
        for(int i=cutCy;i<bitlen;++i){
            cy2[i]=parents[0].getBitCy()[i];
            cy1[i]=parents[1].getBitCy()[i];
        }




        for(int i=cutH.size()-1;i>=0;--i){
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

    vector<double> hillclimb(double x,double phi, Fitness fit, double cx, double cy, vector<double> h){
        double stepX=0.000001;
        double stepPhi=0.0000001;
        
        bool top=0;
        double bestfit=0;

        while(!top){

            




            while(fit(x+stepX,phi,true,cx,cy)>fit(x,phi,true,cx,cy,h)) x+=stepX;
            while(fit(x-stepX,phi,true,cx,cy)>fit(x,phi,true,cx,cy,h)) x-=stepX;
            
            while(fit(x,phi,true,cx+stepX,cy)>fit(x,phi,true,cx,cy,h)) cx+=stepX;
            while(fit(x,phi,true,cx-stepX,cy)>fit(x,phi,true,cx,cy,h)) cx-=stepX;
            while(fit(x,phi,true,cx,cy+stepX)>fit(x,phi,true,cx,cy,h)) cy+=stepX;
            while(fit(x,phi,true,cx,cy-stepX)>fit(x,phi,true,cx,cy,h)) cy-=stepX;

            while(fit(x,phi+stepPhi,true,cx,cy)>fit(x,phi,true,cx,cy,h)) phi+=stepPhi;
            while(fit(x,phi-stepPhi,true,cx,cy)>fit(x,phi,true,cx,cy,h)) phi-=stepPhi;

            if(bestfit-fit(x,phi,true,cx,cy)<0.00000001) top=1;
            bestfit=fit(x,phi,true);
        }

        return {x,phi,cx,cy};
    }
};




namespace plt=matplotlibcpp;

void plotCell(double x, double phi,string plotname="crystall.png",double fitness=0, double cx=0, double cy=0){
    double a=1./sqrt(x*sin(phi)*density);
    vector<double> x1={a,0,0};
    vector<double> x2={a*x*cos(phi),a*x*sin(phi),0};
    vector<double> c={cx,cy};

    vector<double> p1,p2;
    vector<double> u1,u2;
    vector<double> up1,up2;

    
    int lmax=2;
    lmax++;
    int kmax=2;
    kmax++;


    plt::clf();

    for(int i=-layers/2;i<=(layers+1)/2;i++){
    for(int k=-kmax;k<=kmax;++k){
        for(int l=-lmax;l<=lmax;++l){

            switch (layers){
                case 3: {
                    vector<double> upper_vec=sumvec(scalevec(x1,k),scalevec(x2,l),scalevec(c,-1));
                    up1.push_back(upper_vec[0]);
                    up2.push_back(upper_vec[1]);
                }
                case 2: {
                    vector<double> under_vec=sumvec(scalevec(x1,k),scalevec(x2,l),c);
                    u1.push_back(under_vec[0]);
                    u2.push_back(under_vec[1]);
                }
                case 1: {
                    vector<double> vec=sumvec(scalevec(x1,k),scalevec(x2,l));
                    p1.push_back(vec[0]);
                    p2.push_back(vec[1]);
                }
                default: break;
            }
        }
    }
    }
   
    bool bo=plt::plot(p1,p2,"ro",u1,u2,"bs",up1,up2,"gv");
    string title="x="+to_string(x)+" , phi="+to_string(phi/PI*180.)+" , fitness="+to_string(fitness);
    title="d="+to_string(d)+" lambda, with "+to_string(layers)+" layers and density="+to_string(density);
    plt::title(title);
    //cout<<bo;
    //bo=plt::plot(p1,p2);
    //cout<<bo;
    plt::axis("equal");
    //plt::show();
    plt::save(plotname);
    
    cout<<"Plot Saved"<<endl;   
}




int main(){
    int ind_size=1000;

    Fitness fit;
    Generation gen(ind_size,8,6);

  
   
    for(int j=1;j<=10;j++){
        /*d=0.1*j*lambda;
        tric.setH();*/
        for(int i=0;i<100;i++){
            gen_i=i+1;
            //d=j*0.1;
                //high_resolution_clock::time_point t_start_parallel = high_resolution_clock::now();
           // cout<<endl<<"Generation "<<i<<endl;
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
        //gen.printIndiv();
        cout<<endl<<endl<<"Winners of the Evolution Contest "<<j<<" :"<<endl;
        Individual best=gen.getLastBest();
        //best.correctCell();
        cout<<best.lol<<" "<<best.tol<<endl;
        //best.printStats();
        //cout<<endl<<endl<<"ALL THE INDIVS"<<endl;
        //gen.printIndiv();
        Climber climb;
        vector<double> top= climb.hillclimb(best.getX(),best.getPhi(),fit,best.getCx(),best.getCy(),best.getH());
        cout<<endl<<endl<<"After he climbed the hill:"<<endl<<"X="<<top[0]<<endl<<"Phi="<<top[1]*180/PI<<endl<<"a="<<(1./sqrt(density*top[0]*sin(top[1])))<<endl<<"cX="<<top[2]<<endl<<"cY="<<top[3]<<endl<<"Fitness="<<fit(top[0],top[1],false,top[2],top[3],best.getH())<<endl;

        string plotname="crystall";
        plotname+=to_string(j);
        plotname+=".png";
        plotCell(top[0],top[1],plotname,fit(top[0],top[1]),top[2],top[3]);
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