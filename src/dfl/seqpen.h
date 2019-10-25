/* variables, arrays, and routines for seqpen */

  /********************* variables ***********************/
#ifndef SEQPEN_INCLUDED
#define SEQPEN_INCLUDED
class seqpen_optimizer
{
private:
	int nf_max, num_fal, num_iter, num_funct;
	int iprint, istop;
	int i_corr, i_corr_fall;
	double alfa, alfa_max, alfa_stop, obj, violvalue, zobj;

	/********************* arrays **************************/
	/* for dynamic allocation:
	 * (1) uncomment the left hand side definitions, and comment
	 *     out the right hand side definitions
	 * (2) in seqpmain.c: uncomment 'malloc' statements at beginning
	 *     of file and 'free' statements at end of file
	 *
	 * caution: all arrays are used starting with index = 1
	 * hence all dimensions must be increased by 1
	 *
	 *  dynamic alloc     min dimension      static allocation
	 *-------------------------------------------------------*/
	double *alfa_d;
	double *dir;
	double *fstop;

	double *lb;
	double *ub;
	double *xreal;
	double *zvec;

	double *eps;
	double *epsiniz;
	double *constr;

	int nreal; /* number of variables */
	int ncon;  /* number of constraints */

	/*eject*/
	/************************ routines ************************/

	/* seqpenmain.c */
	/* main program */

	/* seqpen.c */
	void seqpen();
	void seqpenmain();
	void displaysolution(char *label);
	double funct(double *x);
	void linesearchbox();
	void teststop();

	/* problem.c */
	void fconstr(double *x);
	double fobj(double *x);
	void setbounds();
	void setdim();
	void startp(double *x);
};

/************* last record of seqpen.h **********/
#endif
