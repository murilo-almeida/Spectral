//
//  Header.h
//  SDG
//
//  Created by Murilo Almeida on 10/26/13.
//  Copyright (c) 2013 Murilo Almeida. All rights resized.
//

#ifndef Tstruct_Header_h
#define Tstruct_Header_h

struct Tstruct {
	int size;
	int * i;
	int * j;
	double * x;
	
	Tstruct() {size =0;};
	
	Tstruct(int n) {
		set(n);
	};
	
	~Tstruct() {
		//cout << "Destruindo Tstruct\n";
		if(size>0) {
			delete [] i; i=nullptr;
			delete [] j; j=nullptr;
			delete [] x; x=nullptr;
		}
	};
	
	void set(int n){
		size=n;
		i = new int [n];
		j = new int [n];
		x = new double [n];
	};
};

struct Vstruct {
	
	int size;
	double * v;
	
	Vstruct() {size = 0;};
	
	Vstruct(int n) { set(n); };
	
	void set(int n) {
		size = n;
		v = new double [size];
	};
	
	~Vstruct() {
		//cout << "Destruindo Vstruct\n";
		if(size>0) {
			delete [] v; v=nullptr;
		}
	};
};


#endif
