#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <initializer_list>
#include <limits>
#include <vector>
#include <tuple>
#include <memory>
#include <queue>
#include <cmath>

using std::queue;
using std::string;
using std::cout;
using std::cerr;
using std::endl;
using std::ofstream;
using std::vector;
using std::tuple;
using std::tie;
using std::make_tuple;
using std::shared_ptr;
#include "io.h"
#include "matrix.h"

#include "MyObject.h"

/*void BFS(Image& res, uint **arr, uint , uint , const uint , const uint , uint);
void FindCentr(Image &, uint** , uint ,const uint , const uint, std::tuple<int, int, bool>*, bool);
int findmin(uint **,const uint ,const uint , const uint , const uint );
int findmax(uint **,uint ,const uint , const uint , const uint );
int num_cogs(uint **, uint , uint , double);
void to_binary(Image &, uint **);
void init_arr(const Image& , uint**);
int comp_BFS(Image &res, uint** arr);
int find_idx( std::tuple<int, int, bool>* , string s[3], uint, Image &, uint**);
void dist_btw_centr(std::tuple<int, int, bool>* ,uint, uint ,uint*,const uint);
void paint(uint , uint , uint , uint , Image &,uint **, uint , uint );
*/
tuple<int, vector<shared_ptr<IObject>>, Image>

//repair_mechanism(const Image& in, string str[3])
repair_mechanism(Image& in, string str[3])
{
	/*uint **arr = new uint*[in.n_rows];
	for(uint i = 0; i < in.n_rows; i++) arr[i] = new uint[in.n_cols];	
	init_arr(in, arr);//arr<- 0|1;
	Image res(in.n_rows, in.n_cols);//create photo;
	to_binary(res, arr);//black&white photo;
	uint L = comp_BFS(res, arr);//count of compononets && BFS
	
	std::tuple<int, int, bool>* card = new std::tuple<int, int, bool>[L];//coordinaty centra
	FindCentr(res, arr, L, res.n_rows, res.n_cols, card, false);
	int result_idx = find_idx(card, str, L, res, arr); //we're HERE
	
	auto object_array = vector<shared_ptr<IObject>>();
	for(uint i=0; i < L; i++){
		uint x, y;
		bool f;
		tie(x, y, f) = card[i];
		if(f){
			object_array.push_back(shared_ptr<IObject>(new Axis(make_tuple(y, x))));
		}
		else{
			int min_r = findmin(arr, x, y, res.n_rows, res.n_cols);
			int max_R = findmax(arr, x, y, res.n_rows, res.n_cols);
			object_array.push_back( shared_ptr<IObject>(new Gear(make_tuple(y, x), min_r, max_R, false, num_cogs(arr, x, y, (min_r+max_R+5)/2))));
		}
	}
	delete [] card;
	delete [] arr;*/
	#define PI 3.3.14159265359;
	for(int i=0; i<3; i++) cout << str[i] << ' ';
	int result_idx = 0;
	auto object_array = vector<shared_ptr<IObject>>();
	for(uint i = 1; i+1 < in.n_rows ; i++){
   		for(uint j=1; j+1 < in.n_cols ; j++){
   			float Gx = 0, Gy = 0, G;
   			int r, g, b;
			tie(r, g, b) = in(i-1, j-1);
			float f1 = 0.299*r + 0.587*g + 0.114*b;
		    tie(r, g, b) = in(i-1, j);
			float f2 = 0.299*r + 0.587*g + 0.114*b;
		    tie(r, g, b) = in(i-1, j+1);
			float f3 = 0.299*r + 0.587*g + 0.114*b;
		    tie(r, g, b) = in(i, j-1);
			float f4 = 0.299*r + 0.587*g + 0.114*b;
		    tie(r, g, b) = in(i, j+1);
			float f6 = 0.299*r + 0.587*g + 0.114*b;
		    tie(r, g, b) = in(i+1, j-1);
			float f7 = 0.299*r + 0.587*g + 0.114*b;
		    tie(r, g, b) = in(i+1, j);
			float f8 = 0.299*r + 0.587*g + 0.114*b;
			tie(r, g, b) = in(i+1, j+1);
      		float f9 = 0.299*r + 0.587*g + 0.114*b;
			Gx = f9 + f7 + 2*f8 - (f1 + 2*f2 + f3);
            Gy = f9 + 2*f6 + f3 - (f1 + 2*f4 + f7);
            G = sqrt(Gx*Gx + Gy*Gy);//modul of gradient
           	if(G > 255) G = 0;
           	if(G < 0) G = 255;
           	in(i, j) = make_tuple(G, G, G);
      	}   			
   }
//    return make_tuple(result_idx, object_array, res.deep_copy());
    return make_tuple(result_idx, object_array, in.deep_copy());
}

int main(int argc, char **argv)
{
    if (argc != 4)
    {
        cout << "Usage: " << endl << argv[0]
             << " <in_image.bmp> <out_image.bmp> <out_result.txt>"<< endl;
        return 0;
    }
	try {
        Image src_image = load_image(argv[1]);
		string s(argv[1]), str[3];
		string s1[3] = {"_1.bmp", "_2.bmp", "_3.bmp"};
		uint i=0; 
		while(s[i] != '.') i++; 
		s.erase(i, 4);
		for(i=0; i < 3; i++){ str[i] = s + s1[i];}
		ofstream fout(argv[3]);
		vector<shared_ptr<IObject>> object_array;
        Image dst_image;
        int result_idx;
        tie(result_idx, object_array, dst_image) = repair_mechanism(src_image, str);
        save_image(dst_image, argv[2]);
        fout << result_idx << endl;
        fout << object_array.size() << endl;
        for (const auto &obj : object_array)
            obj->Write(fout);
    } catch (const string &s) {
        cerr << "Error: " << s << endl;
        return 1;
    }
}/*
int find_idx( std::tuple<int, int, bool>* card, string s[3], uint L,Image & res,uint ** arr1){
	uint x, y;
	bool f;
	for(uint i=0; i < L; i++){
		tie(x, y, f) = card[i];
		if(f){
			for(uint j=0; j < 3; j++){
				char *str = new char[s[j].length() + 1];
				strcpy(str, s[j].c_str());
				Image new_image = load_image(str);
				uint** arr = new uint*[new_image.n_rows];
				for(uint k = 0; k < new_image.n_rows; k++) arr[k] = new uint[new_image.n_cols];	
				init_arr(new_image, arr);//arr<- 0|1;
				to_binary(new_image, arr);//black&white photo;
				std::tuple<int, int, bool>* XY = new tuple<int, int, bool>[1];
				FindCentr(new_image, arr, 1, new_image.n_rows, new_image.n_cols, XY, true);//v XY centr new_image
				uint X, Y;
				bool F;
				tie(X, Y, F) = XY[0];
				int min_r = findmin(arr, X, Y, new_image.n_rows, new_image.n_cols);
				int max_R = findmax(arr, X, Y, new_image.n_rows, new_image.n_cols);//r, R _1|2|3 new_image
				uint * dist = new uint[L-1];
				//cout << x << ' '  << y << " VS " << X << ' ' << Y << endl;
				dist_btw_centr(card, x, y, dist, L);
				uint cnt = 0;
		//		cout << min_r << ' ' << max_R << endl;
				for(uint l=0, k=0; l < L; l++){
					uint x1, y1;
					bool f1;
					tie(x1, y1, f1) = card[l];
					if(!f1){
						uint r = findmin(arr1, x1, y1, res.n_rows, res.n_cols);
						uint R = findmax(arr1, x1, y1, res.n_rows, res.n_cols);
	//					cout << x1 << ':' << y1 << "->" << R+max_R << ':' << r+min_r << ':' << r+max_R << ':' << R+min_r << " = " << dist[k] << endl;
						if((dist[k] < R+max_R) && (dist[k] > r+min_r) && (r+max_R < dist[k]) && (R+min_r < dist[k])) cnt++;
						k++;
						if(cnt == 2) {
				//			cout << j+1 << endl;
							paint(x, y, X, Y, res, arr, new_image.n_cols, new_image.n_rows);//arr1 big photo, arr new kartina
						//	delete [] arr;
							//delete [] str;
							return j+1;
						}
					}
				}
				//sopostavit' i return j;
		//		delete [] arr;
			//	delete [] str;
			}
		}		
	}
	return 0;
}
void paint(uint x, uint y, uint X, uint Y, Image &res,uint** arr, uint n, uint m){
	uint xr = x - X;
	uint yr = y - Y;
	for(uint i=0; i < n; i++){
		for(uint j=0; j < m; j++){
			if(arr[i][j] == 1){
				res(i+xr, j+yr) = make_tuple(0, 0 , 255); 
			}
		}
	}
}

void dist_btw_centr(std::tuple<int, int, bool>* card,uint x0, uint y0,uint* dist,uint L){
	uint x1, y1;
	bool f;
	for(uint i=0, k=0; i < L; i++){
		tie(x1, y1, f) = card[i];
		if(!f){
			dist[k] = round( sqrt((x0-x1)*(x0-x1) + (y0-y1)*(y0-y1)) );
			k++;
		}
	}
}
int comp_BFS(Image & res,uint** arr){
	uint L = 2;
	for (uint j = 0; j < res.n_cols; j++) 
		for (uint i = 0; i < res.n_rows; i++)
	 			if(arr[i][j]==1){BFS(res, arr, i, j,res.n_rows, res.n_cols, L++);}
	return L-2;
}
void init_arr(const Image& in, uint ** arr){
	uint r, g, b;
	for(uint i=0; i < in.n_rows; i++){
		for(uint j=0; j < in.n_cols; j++){
			tie(r, g, b) = in(i, j);
			if(r >= 155 || g >= 155){ arr[i][j] = 1;}
			else {arr[i][j] = 0;}
		}
	}
}
void to_binary(Image & in, uint ** arr){
	for(uint j=0; j < in.n_cols; j++){
		for(uint i=0; i < in.n_rows; i++){
			if(arr[i][j] == 1){in(i, j) = make_tuple(255, 255, 255);}
			else {in(i, j) = make_tuple(0, 0, 0);}
		}
	}
}
int num_cogs(uint **arr,const uint x0,const	uint y0, double r){
	double n = 0, pi = 3.14159265358979323846;
	int cnt = 0;
	bool f = true, t = false;
	int x1 = x0 + r*cos(n);
	int y1 =  y0 + r*sin(n);
	if(arr[x1][y1] == arr[x0][y0])
		if(!(arr[x1+1][y1]+arr[x1-1][y1]+arr[x1][y1+1]+arr[x1][y1-1] < 4*arr[x0][y0])) t = true;
	while (n < 2*pi){
		x1 = x0 + r*cos(n);
		y1 =  y0 + r*sin(n);
		n += 1/r;
		if((arr[x1][y1] == arr[x0][y0])){
	//		res(y1, x1) = make_tuple(0, 0, 255);
			if (f) cnt++;			
			f = false;
		}
		if((arr[x1][y1] != arr[x0][y0])){
//			res(y1, x1) = make_tuple(150, 150, 15);
			f = true;
		}	
	}
	if (t) cnt--;
	return cnt;
}
int findmin(uint **arr,const uint x0,const uint y0, const uint n, const uint m){
	uint x1=0, y1=0;
	for(uint i=0; i < n; i ++){
		for(uint j=0; j < m; j++){
			if((arr[i][j] == arr[x0][y0]) && (i < n-1) && (i > 0)&& (j < m-1) && (j > 0)){//prinadlejit etomu mnojestvu
				if(arr[i+1][j]+arr[i-1][j]+arr[i][j+1]+arr[i][j-1] < 4*arr[x0][y0]){//granichnaya to4ka
					if(sqrt((i-x0)*(i-x0)+(j-y0)*(j-y0)) < sqrt((x1-x0)*(x1-x0)+(y1-y0)*(y1-y0))){//rasstoyanie
						x1 = i;
						y1 = j;
					//	res(i, j) = make_tuple(255, 255, 0);
					}
				}
			}
		}
	}
	return sqrt((x1-x0)*(x1-x0)+(y1-y0)*(y1-y0));
}
int findmax(uint **arr,const uint x0,const uint y0, const uint n, const uint m){
	uint x1=x0, y1=y0;
	for(uint i=0; i < n; i ++){
		for(uint j=0; j < m; j++){
			if((arr[i][j] == arr[x0][y0]) && (i < n-1) && (i > 0)&& (j < m-1) && (j > 0)){//prinadlejit etomu mnojestvu
				if(arr[i+1][j]+arr[i-1][j]+arr[i][j+1]+arr[i][j-1] < 4*arr[x0][y0]){//granichnaya to4ka
					if(sqrt((x1-x0)*(x1-x0)+(y1-y0)*(y1-y0)) < sqrt((i-x0)*(i-x0)+(j-y0)*(j-y0))){//rasstoyanie
						x1 = i;
						y1 = j;
					}
				}
			}
		}
	}
	return sqrt((x1-x0)*(x1-x0)+(y1-y0)*(y1-y0));
}
void FindCentr(Image & res, uint **arr, uint N,const uint n, const uint m, std::tuple<int, int, bool>* card, bool f){
	uint* mas = new uint[N];
	uint* X = new uint[N];
	uint* Y = new uint[N];
	for(uint i=0; i < N; i++){
		X[i] = 0;
		Y[i] = 0;
 		mas[i] = 0;
	}
	uint k;
	if(!f) k = 2; 
	else k = 1;
	for(uint i=0; i < n; i++){
		for(uint j=0; j < m; j++){
			if(arr[i][j] != 0){
				mas[arr[i][j]-k]++;
				X[arr[i][j]-k] += i;
				Y[arr[i][j]-k] += j;
			}
		}
	}
	uint min = mas[0];
	f = false;//Exist Axis???
	for(uint i=0; i < N; i++){
		X[i] /= mas[i];
		Y[i] /= mas[i];
		res(X[i], Y[i]) = make_tuple(0, 0, 0);//vydelenie centra
		if(mas[i] < min){ min = mas[i];}
	}
	for(uint i=0; i < N; i++){
		if((min != mas[i]) && (2 * min < mas[i])){ f = true; break;}
	}
	
	for(uint i=0; i < N; i++){
		if((mas[i] == min) && f){card[i] = std::make_tuple(X[i], Y[i], true);}
		else{card[i] = std::make_tuple(X[i], Y[i], false);}
	}
}
void BFS(Image& res, uint **arr, uint x, uint y, const uint N, const uint M, uint L){
	queue<uint> queueX;
	queue<uint> queueY;
	queueX.push(x);
	queueY.push(y);
	while (!queueX.empty()){
		y = queueY.front();
		x = queueX.front();
		if ((x+1 < N) && (arr[x+1][y]==1)){ queueX.push(x+1); queueY.push(y); arr[x+1][y] = L; res(x+1, y) = make_tuple((300*L+100)%256, (300*L+100)%256, (300*L+100)%256);};
		if ((y+1 < M) && (arr[x][y+1]==1)){ queueX.push(x); queueY.push(y+1); arr[x][y+1] = L; res(x, y+1) = make_tuple((300*L+100)%256, (300*L+100)%256, (300*L+100)%256);};
		if ((x > 0 && x < N) && (arr[x-1][y]==1)){queueX.push(x-1); queueY.push(y); arr[x-1][y] = L; res(x-1, y) = make_tuple((300*L+100)%256, (300*L+100)%256, (300*L+100)%256);};
		if ((y > 0 && y < M) && (arr[x][y-1]==1)){queueX.push(x); queueY.push(y-1); arr[x][y-1] = L; res(x, y-1) = make_tuple((300*L+100)%256, (300*L+100)%256, (300*L+100)%256);};
		queueY.pop();
		queueX.pop();
	}
}*/