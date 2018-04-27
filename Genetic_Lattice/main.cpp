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

//double const a=1;
double const lambda=1;
double const d=lambda*1;
int n=100;

int const bitlen=8;

double rcut=lambda*13;
//int kmax=15;
//int lmax=15;

int parent_amount=2;        //changes in code have to be made!!!
double mutation_rate=0.01;

vector<double> veczero={0,0,0};

int threadcount=8;

struct Tric{
    public:
    double const phi=PI/3.;
    double const x=1;
    double const cx=0;
    double const cy=0;
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
    for(int i=0; i<x1.size();++i) x.push_back(x1.at(i)+x2.at(i)+x3.at(i));
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
    double operator() (double x, double phi, Potential v,bool min=false, double cx=0, double cy=0) const{
        double sum=0;
        double a=1./sqrt(x*sin(phi));
        vector<double> x1={a,0,0};
        vector<double> x2={a*x*cos(phi),a*x*sin(phi),0};
        vector<double> x3={cx,cy,d};

        if(min){
        vector<vector<double>> cell=minimizeCell(x1,x2);   
        x1=cell[0];
        x2=cell[1];
        }
        
        int lmax=rcut/a/x;
        lmax++;
        int kmax=rcut/a;
        kmax++;

        

        for(int k=-kmax;k<=kmax;++k){
            for(int l=-lmax;l<=lmax;++l){
                if(l!=0||k!=0){
                    sum+=v(euclid(sumvec(scalevec(x1,k),scalevec(x2,l),veczero)));            //optimizable
                    sum+=v(euclid(sumvec(scalevec(x1,k),scalevec(x2,l),x3)));
                    //if(euclid(sumvec(scalevec(x1,k),scalevec(x2,l)))==0) cout<<"Null found at kl"<<k<<" "<<l<<endl;
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
            tric_sum=latticeSum(tric.x,tric.phi,potential,tric.cx,tric.cy);
        }
        
        double operator() (double x,double phi,bool min=false,double cx=0,double cy=0) const{
            //if(phi<PI/6) return 0;
            double sum=latticeSum(x,phi,potential,min,cx,cy);
            return exp(1-sum/tric_sum);
        }
} fit;





class Individual{
    vector<bool> x;
    vector<bool> phi;
    bitset<bitlen> cx;
    bitset<bitlen> cy;

    double fitness=0;

    public:

    vector<bool> getVecX(){
        return x;
    }
    vector<bool> getVecPhi(){
        return phi;
    }

    bitset<bitlen> getBitCx(){
        return cx;
    }

    bitset<bitlen> getBitCy(){
        return cy;
    }

    void setVecX(double dx){
        dx*=pow(2,x.size());
        int ix=--dx;
        vector<bool> vecX;

        for(int i=0;i<x.size();i++){
            if(ix){
                vecX.push_back(ix&1);
                ix>>=1;
            }
            else vecX.push_back(0); 
        }
        
        //reverse(vecX.begin(),vecX.end());
        x=vecX;
    }

    void setVecPhi(double dx){
        dx*=pow(2,phi.size());
        dx/=PI/2;
        int ix=--dx;
        
        vector<bool> vecX;

        for(int i=0;i<phi.size();i++){
            if(ix){
                vecX.push_back(ix&1);
                ix>>=1;
            }
            else vecX.push_back(0); 
        }
        
        //reverse(vecX.begin(),vecX.end());
        phi=vecX;
    }

    void setBitCx(double dx){
        dx*=pow(2,cx.size());
        int ix=--dx;

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
        int ix=--dx;

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



    double getX(){
        double accum=accumulate(x.rbegin(),x.rend(),0,[](int i, int j){ return (i<<1)+j;});    //shift the x bits to the left, so that y can be inserted on right site
        accum++;
        accum/=pow(2,x.size());
        
        return accum;
    }
    double getPhi(){
        double accum=accumulate(phi.rbegin(),phi.rend(),0,[](int x, int y){ return (x<<1)+y;});
        accum++;
        accum*=PI/2;
        accum/=pow(2,phi.size());
        return accum;
    }


    //not sure if right
    double getCx(){
        double accum=cx.to_ulong();
        accum++;
        accum*2;
        accum/=pow(2,cx.size());
        return accum;
    }
    
    double getCy(){
        double accum=cy.to_ulong();
        accum++;
        accum*2;        //to make max length to 2
        accum/=pow(2,cy.size());
        return accum;
    }


    double getFitness(){  
        return fitness;
    }

    void calcFitness(){
        fitness=fit(getX(),getPhi(),false,getCx(),getCy());
    }

    void printStats(){
        cout<<"Fitness="<<fitness<<endl<<"x="<<getX()<<endl<<"phi="<<getPhi()/PI*180<<endl<<"a="<<1./sqrt(sin(getPhi())*getX())<<endl<<"cx="<<getCx()<<endl<<"cy="<<getCy()<<endl;
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
        for(int i=0;i<x_length;i++) x.push_back(dis(gen));
        for(int i=0;i<phi_length;i++) phi.push_back(dis(gen));
        for(int i=0;i<cx.size();i++) cx[i]=dis(gen);
        for(int i=0;i<cy.size();i++) cy[i]=dis(gen);
    }
    //child Creator
    Individual(vector<bool> xVec, vector<bool> phiVec, bitset<bitlen> cx, bitset<bitlen> cy){ //has to be changed for more or less than 2 parents!!!!!!
        x=xVec;
        phi=phiVec;
        this->cx=cx;
        this->cy=cy;
        mutate();
        correctCell();
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

    }


    void correctCell(){
        double x=getX();
        double phi=getPhi();

        double a=1./sqrt(x*sin(phi));
        vector<double> x1={a,0};
        vector<double> x2={a*x*cos(phi),a*x*sin(phi)};
        
        
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
        x=euclid(x2)/a;
        phi=acos(x2[0]/a/x); 

        setVecX(x);
        setVecPhi(phi);
    }




    vector<vector<double>> minimizeCell(vector<double> x1, vector<double> x2) const{
        bool minimal=0;
        double u=circ(x1,x2);
        while(!minimal){
            bool xchange=0;
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
            t.push_back(thread(&Generation::createChildren,this,&children,size/threadcount));       //size has to be devisable by threadcoutn
        }
        for(int i=0;i<threadcount;++i){
            t.at(i).join();
        }
        
        //cout<<"Children become Parents"<<endl;
        indiv=children;
        children.clear();
        //printBest();
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
        
        vector<bool> x1(parents.at(0).getVecX().size());
        vector<bool>phi1(parents.at(0).getVecPhi().size());
        vector<bool>x2(parents.at(0).getVecX().size());
        vector<bool>phi2(parents.at(0).getVecPhi().size());
        
        bitset<bitlen> cx1;
        bitset<bitlen> cx2;
        bitset<bitlen> cy1;
        bitset<bitlen> cy2;
        
        /*using namespace std::chrono;
        high_resolution_clock::time_point t_start_parallel = high_resolution_clock::now();
        *///cout<<"Genes will be created!"<<endl;
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
        
        /*high_resolution_clock::time_point t_end_parallel = high_resolution_clock::now();
        duration<double> time_parallel = t_end_parallel - t_start_parallel;
        cout << "Execution time: " << time_parallel.count()<<endl;
        */

        //cout<<"The two children will be created"<<endl;
        Individual child1(x1,phi1,cx1,cy1);
        Individual child2(x2,phi2,cx2,cy2);

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

    vector<double> hillclimb(double x,double phi, Fitness fit, double cx, double cy){
        double stepX=0.001;
        double stepPhi=0.00001;
        
        bool top=0;
        double bestfit=0;

        while(!top){
            while(fit(x+stepX,phi,true,cx,cy)>fit(x,phi,true,cx,cy)) x+=stepX;
            while(fit(x-stepX,phi,true,cx,cy)>fit(x,phi,true,cx,cy)) x-=stepX;
            while(fit(x,phi+stepPhi,true,cx,cy)>fit(x,phi,true,cx,cy)) phi+=stepPhi;
            while(fit(x,phi-stepPhi,true,cx,cy)>fit(x,phi,true,cx,cy)) phi-=stepPhi;
            while(fit(x,phi,true,cx+stepX,cy)>fit(x,phi,true,cx,cy)) cx+=stepX;
            while(fit(x,phi,true,cx-stepX,cy)>fit(x,phi,true,cx,cy)) cx-=stepX;
            while(fit(x,phi,true,cx,cy+stepX)>fit(x,phi,true,cx,cy)) cy+=stepX;
            while(fit(x,phi,true,cx,cy-stepX)>fit(x,phi,true,cx,cy)) cy-=stepX;
            if(bestfit-fit(x,phi)<0.01) top=1;
            bestfit=fit(x,phi,true);
        }

        return {x,phi,cx,cy};
    }
};




namespace plt=matplotlibcpp;

void plotCell(double x, double phi,string plotname="crystall.png",double fitness=0, double cx=0, double cy=0){
    double a=1./sqrt(x*sin(phi));
    vector<double> x1={a,0};
    vector<double> x2={a*x*cos(phi),a*x*sin(phi)};
    vector<double> c={cx,cy};

    vector<double> p1,p2;
    vector<double> u1,u2;

    
    int kmax=5;
    int lmax=5;


    for(int k=-kmax;k<=kmax;++k){
        for(int l=-lmax;l<=lmax;++l){
            vector<double> vec=sumvec(scalevec(x1,k),scalevec(x2,l));
            vector<double> under_vec=sumvec(scalevec(x1,k),scalevec(x2,l),c);
            p1.push_back(vec[0]);
            p2.push_back(vec[1]);
            u1.push_back(under_vec[0]);
            u2.push_back(under_vec[1]);
        }
    }
    plt::clf();
    bool bo=plt::plot(p1,p2,"ro",u1,u2,"bs");
    string title="x="+to_string(x)+" , phi="+to_string(phi/PI*180.)+" , fitness="+to_string(fitness);
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
    Generation gen(ind_size,8,8);

  
    
    for(int j=0;j<10;j++){
        for(int i=0;i<100;i++){
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
        gen.printIndiv();
        cout<<endl<<endl<<"Winner of the Evolution Contest "<<j<<" :"<<endl;
        Individual best=gen.getLastBest();
        Climber climb;
        vector<double> top= climb.hillclimb(best.getX(),best.getPhi(),fit,best.getCx(),best.getCy());
        cout<<endl<<endl<<"After he climbed the hill:"<<endl<<"X="<<top[0]<<endl<<"Phi="<<top[1]*180/PI<<endl<<"cX="<<top[2]<<endl<<"cY="<<top[3]<<endl<<"Fitness="<<fit(top[0],top[1])<<endl;

        string plotname="crystall";
        plotname+=to_string(j);
        plotname+=".png";
        plotCell(top[0],top[1],plotname,fit(top[0],top[1]),top[2],top[3]);
    }
    
    
    cout<<"FINISHED! \nEnter anything to close"<<endl;
    bool b;
    cin>>b;
    return 0;
}