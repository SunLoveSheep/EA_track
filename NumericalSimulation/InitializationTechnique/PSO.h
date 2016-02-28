//PSO parameters and operators
#ifndef _PSO_H
#define _PSO_H

struct PSOparameter
{
	double c1,c2,w;//updating coefficients
	//double *pbest;//to store the best result ever found by each particle
	double gbest;//global best for the iteration
};

class PSO
{
public:
	PSO();
	virtual~PSO();

	double *pbest;
	double **pbestposition;//the pbest position
	double *gbestposition;//position information for global best
	double **v;//velocity matrix
	double *vmin,*vmax;//minimum and maximum of velocity

	int FE;

	void PSOparticleupdate();
	void PSOpbestupdate();//after every iteration, update pbest
	void PSOgbestupdate();
	void Initial();

	void PSOprocess(int loop);
};

#endif