#ifndef RAYSOLVER_H
#define RAYSOLVER_H


#include <vector>
#include <iostream>
#include "RayTrace.h"
#include "RayTrace_IceModels.h"
#include "Vector.h"
#include <fstream>
#include <iomanip>
#include <boost/lexical_cast.hpp>
#include <boost/program_options.hpp>

class Position;
class IceModel;

//--------------------------------------------------
// struct underline{
// 	static char esc;
// 	const std::string& str;
// 	underline(const std::string& s):str(s){}
// 	friend std::ostream& operator<<(std::ostream& os, const underline& u);
// };
// 
// char underline::esc=0x1B;
// 
// std::ostream& operator<<(std::ostream& os, const underline& u){
// 	return(os << underline::esc << "[4m" << u.str << underline::esc << "[0m");
// }
//-------------------------------------------------- 

template<typename positionType>

class pathPrinter{
public:
	pathPrinter(){}
	void operator()(const positionType& p, RayTrace::RKStepType stepType){
		std::cout << p.x << ' ' << p.z << '\n';
	}
};

class RaySolver {

    private:
        void Earth_to_Flat_same_depth(Position &source, Position &target, IceModel* antarctica);
        void Earth_to_Flat_same_angle(Position &source, Position &target, IceModel* antarctica);

    public:
        RaySolver();
//        RaySolver(int argc, char* argv[]);

        void test1();
        void Solve_Ray_org(Position &source, Position &target, std::vector < std::vector <double> > &outputs);
        void Solve_Ray(Position &source, Position &target, IceModel *antarctica, std::vector < std::vector <double> > &outputs);

        int source_over_surface;
        int solution_toggle;    // no solution : 0  solution exist : 1
};


#endif //RAYSOLVER_H
