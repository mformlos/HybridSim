/*
 * Analysis.h
 *
 *  Created on: Apr 28, 2015
 *      Author: maud
 */

#ifndef LIB_ANALYSIS_H_
#define LIB_ANALYSIS_H_

#include <ostream>

template<class ... Args>
class Analysis {
public:
	virtual ~Analysis() {}
	virtual void operator() (const Args& ... args){}
	virtual double value(){return -1.0;}
	virtual std::ostream& print_result(std::ostream& os) {return os;}
};


#endif /* LIB_ANALYSIS_H_ */
