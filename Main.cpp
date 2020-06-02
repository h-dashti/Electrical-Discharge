//	We develop a model of an electrical discharge in inhomogenous isulator. (for example lighting ). When an elctrical discharging occurs
//	the elctrical potential phi, staistfy Laplace equation( Laplaican(phi) = 0 ).  The model specified by these steps:
//
/*
	1- Consider a large boundary circle of radius R and place charge source at origin. Choose the potential phi=0 at the origin(ocuiped site)
	and phi=1 for sites on circumference of the circle. The radius R should be larger than the radius of the growing pattern. 

	2- We using Relaxation method to compuate phi_i for empty sites in circle.

	3- Assign a random number r to each empty in circle. The random number r_i at site i represents a breakdown coefficent and the random
	inhomegenus nature of insulator.

	4- The perimeter sites are the neigbour sites of the discharged pattern(occuiped sites). We form the product r_i*(phi_i)^a for each perimeter
	site i , where a is adjustable parameters.

	5- The perimeter site with maximum value of r*phi^a , breaks down(means its potential equal to zero ).

	6- We use the Relaxation method to recaclulate the value of potential at the remaining unoccuiped sites, and then we repeat steps (4)-(6).




// Site(x,y) = 1	:	represnet occupied site.
// Site(x,y) = 2	:	represnet a perimeter site.
// Site(x,y) = 0	:	represnet an un tested site(empty).

--------------------------------------------------------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------------------------------------------------------
*/


#include <iostream>
#include <fstream>
#include <ctime>
#include <cmath>
#include <cstdlib>
using namespace std;


// define structures
struct Location
{
    int x, y;
    Location() : x(0), y(0)
    {}
    Location(int X, int Y) :x(X), y(Y)
    {}

    Location& operator = (const Location& loc)
    {
        x = loc.x;
        y = loc.y;
        return *this;
    }
};
// End define Structers

//----------------------------------------------------------------------------------- Define Variables
const int R = 100 , L = 2*R+1, R2 = R*R; // radius of circle
double phi[L][L]; // potentail array 
double phiTemp[L][L];
char site[L][L];  // represent that the a site is emtpy or occuiped or perimeter
Location peremiter[L*L/2]; // the array that represnt perimeter location
Location seeds[L*L];//		// the array that represnt occuiped location

int n_seeds;// number of occuiped sites
int n_per; // number of perimeter sites
const double tolerance = 1E-5;	// tolerance for calculation potential
const int pot_iterations = 500;	// number of pot_iterations for solving potential

int nX[] = {1, 0, -1, 0};  // Set Up direction vectors for latice
int nY[] = {0, 1, 0, -1};

double r_Sites[L][L]; // array that repersent the breakdown coefficient for ecah site

const double a = 1.0/4;
const int n_discharging = 200; // the number of site sthat must been discharged
//-----------------------------------------------------------------------------------End Define Variables

//-----------------------------------------------------------------------------------Prototype Methods
inline int Rand( int nMin, int nMax );
inline double Rand( );


void Initialize();
void SolvePotential();
void CalculatePot(double v1[L][L], double v2[L][L], double& diff);

void DischargeDielectric();
void Discharge(int i_max);

void OutPutData()
{
	int i, j, di, dj;
	int count = 0;
	ofstream out ("phi.txt", ios::out );
	for(i=0; i<L; i++)
	{
		di = i-R;
		for(j=0; j<L; j++)
		{
			dj = j-R;
			if(di*di + dj*dj < R2 )
			{
				if(count %3 == 0)
					out << i << " " << j << " " << phi[i][j] << endl ;
				count++;
			}
		}

	}
	out.close();
};


//-----------------------------------------------------------------------------------



int main()
{
	srand((unsigned)time(NULL));

	Initialize();
	DischargeDielectric();
	cout << "\a\n";
	OutPutData();
	return 0;
}

void Initialize()
{
	int i, j, di, dj ;

	for( i = 0; i < L; i++)
		for(j = 0; j < L; j++)
			site[i][j] = 0;

	


	for(i = 0; i < L; i++) // Initilaize all arrays
	{
		di = i-R;
		for(j = 0; j < L; j++)
		{
			dj = j-R;
			if (di*di + dj*dj >= R2)
				phi[i][j] = phiTemp[i][j] = 1;  // set the boundary and out on, equal 1
			else
				phi[i][j] = phiTemp[i][j] = 0; // set inside boundary equal 0

		}
	}
	
	
	phi[R][R] =  phiTemp[R][R] = 0;  // set the center 0, and add it to occuiped(seed) sites.
	n_seeds = 1;
	seeds[0] = Location(R, R);
	site[R][R] = 1;
	

	n_per = 4; // set the peremiter sites
	int nx, ny;
	for(i = 0; i < 4; i++)
	{
		nx = R+nX[i];
		ny = R+nY[i];
		peremiter[i] = Location(nx, ny); 
		site[nx][ny] = 2;
	}



	// randomize each site and calculate r for each site
	for(i = 0; i < L; i++)
		for ( j =0; j< L; j++)
			r_Sites[i][j] = Rand();


	SolvePotential(); // now we find potential



}

void SolvePotential()
{
	/*for(int i = 0; i < L; i++)
		for(int j = 0; j < L; j++)
			phiTemp[i][j] = phi[i][j];*/
		
	double diff = 1;
	for(int i = 0; i < pot_iterations || diff >= tolerance; i++)
	{
		CalculatePot(phi, phiTemp, diff);
		CalculatePot(phiTemp, phi, diff);
	}
	
}
void CalculatePot(double v1[L][L], double v2[L][L], double& diff)
{
	double value;
	diff = 0;
	int i, j, di, dj ,count = 0;
	int j1, j2;

	for(i = 1; i < L-1; i++ )
	{
		di = i-R;
		j1 = (int)(R - sqrt((double)(R2 - di*di)));
		j2 = (int)(R + sqrt((double)(R2 - di*di)));
		for ( j = j1; j <= j2; j++)
		{
			if(site[i][j] == 1) // we dont need to calculate potential of discharged sites (occuiped site)
				continue; 
			dj = j-R;
			if(di*di + dj*dj >= R2) //we dont need to calculate potential on boundary (occuiped site)
				continue;
			value = v1[i][j];
			v2[i][j] = 0.25 * (v1[i + 1][j] + v1[i - 1][j] + v1[i][j + 1] + v1[i][j - 1]);
			diff += abs(value - v2[i][j] );
			count++;
		}
	}

	diff /= count;

}





void DischargeDielectric()
{
	int i, x, y,  i_max;
	double value, max;
	int temp = n_discharging / 100;
	
	cout << "0%" ;
	for(int j = 1; j <= n_discharging; j++)
	{
		max = -1;
		for(i = 0; i < n_per; i++)
		{
			x = peremiter[i].x;
			y = peremiter[i].y;
			value = r_Sites[x][y] * phi[x][y]; //pow(phi[x][y], a);
			if(value > max)
			{
				max = value;
				i_max = i;
			}
		}
		Discharge(i_max);
		SolvePotential();
		if(j % temp == 0)
			cout << "\b\b\b" << (100 * j) / n_discharging << "%";
	}



}

void Discharge(int i_per)
{
	int x = peremiter[i_per].x;
	int y = peremiter[i_per].y;
	site[x][y] = 1;
	phi[x][y] = phiTemp[x][y] = 0; // discharge the target site
	n_seeds++; // now we increase the number of occuiped site and decrease the perimeter number
	seeds[n_seeds-1].x = x;
	seeds[n_seeds-1].y = y;
	peremiter[i_per] = peremiter[n_per - 1];
	n_per--;
		

	int i, xNew, yNew, dx, dy;
	for(i=0; i<4; i++) // now generate new perimeter sites
	{
		xNew = x + nX[i];
		yNew = y + nY[i];
		dx = xNew - R;
		dy = yNew - R;
		if(dx*dx + dy*dy >= R2)  
			continue; 
		if(site[xNew][yNew] == 0) // we can add to perimeter list if the site is empty
		{
			n_per++;
			peremiter[n_per-1].x = xNew;
			peremiter[n_per-1].y = yNew;
			site[xNew][yNew] = 2;
		}
	}
}




inline int Rand( int nMin, int nMax )
{
    int nRange = nMax - nMin;
    int nNum = rand() % nRange;
    return( nNum + nMin );
}
inline double Rand()
{
	return (double)rand() / RAND_MAX;

}