#include "PhononFrequency.h"

#include <gsl/gsl_matrix_complex_double.h>
#include <gsl/gsl_eigen.h>

#include <cmath>
#include <algorithm>

void PhononFrequency(const GeometryVector & k, const Configuration & crystal, PairPotential & pot, std::vector<double> & result2, double Mass)
{
	int dim=crystal.GetDimension();
	int nbr=crystal.NumParticle();
	//Configuration::particle * start = crystal.GetParticle(0);

	gsl_matrix_complex * mat = gsl_matrix_complex_calloc(nbr*dim, nbr*dim);

	for(int i=0; i<nbr; i++)
	{
		AtomInfo a1 = crystal.GetCharacteristics(i);
		const Configuration * pc = &crystal;
		crystal.IterateThroughNeighbors(i, pot.Rcut,
			[&k, &dim, &mat, &i, &pot, &a1, &pc, &nbr] (const GeometryVector & shift, const GeometryVector & LatticeShift, const signed long * PeriodicShift, const size_t srcAtom) ->void
			{
				if(shift.Modulus2()< ::LengthPrecision* ::LengthPrecision)
					return;
				int j=srcAtom;
				AtomInfo a2 = pc->GetCharacteristics(j);
				double ex=std::cos(k.Dot(LatticeShift));
				double iex=std::sin(k.Dot(LatticeShift));
				for(int m=0; m<dim; m++)
				{
					for(int n=0; n<dim; n++)
					{
						if (j < nbr)
						{
							(*gsl_matrix_complex_ptr(mat, i*dim + m, j*dim + n)).dat[0] += ex* pot.SecondDerivative(shift, m, n, a1, a2);
							(*gsl_matrix_complex_ptr(mat, i*dim + m, j*dim + n)).dat[1] += iex* pot.SecondDerivative(shift, m, n, a1, a2);
						}
						(*gsl_matrix_complex_ptr(mat, i*dim + m, i*dim + n)).dat[0] -= pot.SecondDerivative(shift, m, n, a1, a2);
					}
				}
			}
		);
	}

	gsl_eigen_herm_workspace * wor = gsl_eigen_herm_alloc(dim*nbr);
	gsl_vector * eval = gsl_vector_alloc(dim*nbr);

	gsl_eigen_herm(mat, eval, wor);
	result2.clear();
	for(size_t i=0; i<dim*nbr; i++)
		result2.push_back((-1)*gsl_vector_get(eval, i)/Mass);
	std::sort(result2.begin(), result2.end());

	gsl_vector_free(eval);
	gsl_eigen_herm_free(wor);
	gsl_matrix_complex_free(mat);
}

