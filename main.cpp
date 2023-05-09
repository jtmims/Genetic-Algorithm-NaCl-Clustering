#include <iostream>
#include <math.h>
#include <cstdlib>
#include <random>
#include <vector>
#include <algorithm>
#include <fstream>

//CONSTANTS AND GLOBAL VARIABLES
const int GENERATIONS = 100;
const int INITPOPSIZE = 10000;
const double KEEPRATIO = 0.5;
const int POPSIZE = INITPOPSIZE*KEEPRATIO;
const double ALPHA = 0.25;
const double MUTATIONAMOUNT = 20;
const double MUTATIONPERCENT = 0.01;
const int AVGFITCUTOFF = 2000;
int numAtoms;
int numNa;
int numCl;
/*
Initializing random device and various distributions.
Molecule bond lengths should be along the scale of magnitude for an Angstrom making our random guesses fall in the interval (0.01,10).
*/
std::random_device rd;
std::uniform_real_distribution<float> dist(-10,10);
std::uniform_real_distribution<float> perc(0.0,1.0);
std::uniform_real_distribution<float> mut((1-MUTATIONPERCENT),(1+MUTATIONPERCENT));
/*
Common activation function used in Neural Networks.
This function penalizes all negative members the same and has a high slope towads positive members.
This might cause problems for a problem with a large solution space, but this is not an issue for
us as explained above.
The prefactor was chosen arbitrarily.
*/
float expLinUnit(float x){
	if(x>=0){
		return (exp(x)-1);
	}else{
		return x;
	}
}
/*
Function that returns the Euclidean norm in R^3
*/
float euclidNorm(float x_1, float y_1, float z_1, float x_2, float y_2, float z_2){
	float dx = x_1-x_2;
	float dy = y_1-y_2;
	float dz = z_1-z_2;
	return sqrt((dx*dx)+(dy*dy)+(dz*dz));
}
/*
Structure that contains the name of the atom and position.
*/
struct atom{
	std::string name;
	float x,y,z;
};
/*
Structure that will contain all the information for each member of the population.
This includes a vector of atoms and the cooresponding fitness score.
The potential is in terms of eVs and Angstroms for ease of use.
*/
struct member{
	float fitness;
	std::vector<atom> particles;
	//Calculate potentials for each atom and calculate fitness for this member
	void foo(){
		float potential = 0;
		for(int i=0;i<numAtoms;i++){
			for(int j=i+1;j<numAtoms;j++){
				bool kronecker = (particles[i].name != particles[j].name);
				int eta = pow(-1,kronecker);
				float r = euclidNorm(particles[i].x,particles[i].y,particles[i].z,particles[j].x,particles[j].y,particles[j].z);
				if(r <= 0.1){
					potential += 9999;
				}else{
					//std::cout<<r<<std::endl;
					float ep = 55.26349406;
					float k = 1/(4*M_PI*ep*pow(10,-4));
					potential += eta*(k/r) + 1.09*pow(10,3)*exp((-r)/0.321);// + pow((0.1/r),12);
				}
			}
		}
		fitness = expLinUnit(-potential);
	}
};
/*
Function that generates the initial population. 
*/
std::vector<member> generatePopulation(){
	std::vector<member> p;
	for(int i = 0; i < INITPOPSIZE; i++){
		std::vector<atom> a;
		for(int j=0;j<numNa;j++){a.emplace_back(atom{"Na"});}
		for(int j=0;j<numCl;j++){a.emplace_back(atom{"Cl"});}
		//generating random positions with some added rules to decrease the degress of freedom of the system
		for(int j=0;j<numAtoms;j++){
			switch(j){
				case 0:
					a[j].x, a[j].y, a[j].z = 0;
					break;
				case 1:
					a[j].x = dist(rd);
					a[j].y, a[j].z = 0;
					break;
				case 2:
					a[j].x = dist(rd);
					a[j].y = dist(rd);
					a[j].z = 0;
					break;
				default:
					a[j].x = dist(rd);
					a[j].y = dist(rd);
					a[j].z = dist(rd);
			}
		}
		p.emplace_back(member{0,a});
	}
	for(auto& m : p){m.foo();}
	return p;
}
/*
Implementation of the BLX-alpha crossover method of real values.
This function returns a random number in the range [min-alpha*range, max+alpha*range].
*/
float BLXalphaCross(float q_1, float q_2){
	float u, l;
	if(q_1 >= q_2){
		u = q_1;
		l = q_2;
	}else if(q_2 > q_1){
		u = q_2;
		l = q_1;
	}
	float lb = l-ALPHA*(u-l);
	float ub = u+ALPHA*(u-l);
	return (lb + perc(rd)*(ub-lb));
}
/*
Function that handles both selection and crossover.
The roulette wheel selection method is used.
The best member of the population does not mate to change, but does mate with others.
*/
void selectionCrossover(std::vector<member>& p){
	float sum = 0;
	for(int i=0;i<POPSIZE;i++){
		if(p[i].fitness >=0){sum += p[i].fitness;}
	}
	//std::cout<<sum<<std::endl;
	for(int i=1;i<POPSIZE;i++){
		float rand = perc(rd)*sum;
		float check = 0;
		//std::cout << rand << std::endl;
		for(int j=0;j<POPSIZE;j++){
			check += p[j].fitness;
			if(check >= rand){
				for(int a=0;a<numAtoms;a++){
					switch(a){
						case 0:
							p[i].particles[a].x = 0;
							p[i].particles[a].y = 0;
							p[i].particles[a].z = 0;
							break;
						case 1:
							p[i].particles[a].x = BLXalphaCross(p[j].particles[a].x,p[i].particles[a].x);
							p[i].particles[a].y = 0;
							p[i].particles[a].z = 0;
							break;
						case 2:
							p[i].particles[a].x = BLXalphaCross(p[j].particles[a].x,p[i].particles[a].x);
							p[i].particles[a].y = BLXalphaCross(p[j].particles[a].y,p[i].particles[a].y);
							p[i].particles[a].z = 0;
							break;
						default:
							p[i].particles[a].x = BLXalphaCross(p[j].particles[a].x,p[i].particles[a].x);
							p[i].particles[a].y = BLXalphaCross(p[j].particles[a].y,p[i].particles[a].y);
							p[i].particles[a].z = BLXalphaCross(p[j].particles[a].z,p[i].particles[a].z);
					}
				}
				p[i].foo();
				break;
			}
		}
	}
}
/* 
Sort function from the algorithms header file put into a function to keep from writing multiple times.
*/
void sort(std::vector<member>& p){
	std::sort(
    	p.begin(), 
    	p.end(), 
    	[](const auto& l,const auto& r){
     		return l.fitness > r.fitness;
    	});
}
/*
This function mutates a random number of members, besides the best, depending on the constant MUTATIONPERCENTAGE.
*/
void mutate(std::vector<member>& p){
	for(int i=0;i<MUTATIONAMOUNT;i++){
		int rand = 1+perc(rd)*(POPSIZE-1);
		for(int j=0;j<numAtoms;j++){
			switch(j){
				case 0:
					break;
				case 1:
					p[rand].particles[j].x *= mut(rd);
					break;
				case 2:
					p[rand].particles[j].x *= mut(rd);
					p[rand].particles[j].y *= mut(rd);
					break;
				default:
					p[rand].particles[j].x *= mut(rd);
					p[rand].particles[j].y *= mut(rd);
					p[rand].particles[j].z *= mut(rd);
			}
		}
	}
}

/*
Outputs data to console and text files.
*/
void recordData(std::vector<member>& p,int g,std::ofstream &xyz,std::ofstream &consolelog){
	float avg = 0;
	int num = 0;
	xyz << numAtoms << std::endl;
	xyz << "test" << std::endl;
	for(int i=0;i<numAtoms;i++){
		xyz << p[0].particles[i].name << "\t" << p[0].particles[i].x << "\t" << p[0].particles[i].y << "\t" << p[0].particles[i].z << std::endl;
	}
	for(auto& m : p){
		avg += m.fitness;
		num++;
	}
	avg = avg/num;
	std::cout << "Generation: " << g << std::endl;
		std::cout << "Avg Fitness: " << avg << std::endl;
}

int main(){
	//Initializing output files
	std::ofstream xyz;
	xyz.open("config.xyz");
	std::ofstream consolelog; 
	consolelog.open("consolelog.txt");

	//Taking user input to determine the desired system
	std::vector<atom> a;
	std::cout << "Number of Sodium atoms: " << std::endl;
	std::cin >> numNa;
	std::cout << "Number of Chlorine atoms: " << std::endl;
	std::cin >> numCl;
	numAtoms = numNa + numCl;
	
	//Generating a random population, culling the worst members, and outputting the best initial configuration
	std::vector<member> initPopulation = generatePopulation();
	sort(initPopulation);
	std::vector<member> population;
	std::copy(initPopulation.begin(), initPopulation.begin()+POPSIZE, std::back_inserter(population));
	recordData(population, 0, xyz, consolelog);

	for(int g=1;g<=GENERATIONS;g++){
		sort(population);
		selectionCrossover(population);
		//for(auto& m : Population){std::cout << m.fitness << " " << m.r << std::endl;}
		mutate(population);
		sort(population);
		recordData(population,g,xyz,consolelog);
	}

	//std::cout << population[0].fitness << " " << population[0].particles[0].x << " "<< population[0].particles[1].x<< " " << population[0].particles[2].x<< std::endl;
	xyz.close();
	consolelog.close();
}