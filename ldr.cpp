/*
  OM*
  Genome-wide Efficient Mixed Model Association (GEMMA)
Copyright (C) 2011  Xiang Zhou

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
  
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <ctime>
#include <cstring>
#include <algorithm>

#include "gsl/gsl_vector.h"
#include "gsl/gsl_matrix.h"
#include "gsl/gsl_linalg.h"
#include "gsl/gsl_blas.h"
#include "gsl/gsl_eigen.h"
#include "gsl/gsl_randist.h"
#include "gsl/gsl_cdf.h"
#include "gsl/gsl_roots.h"
#include "Eigen/Dense"

#include "gsl/gsl_rng.h"


#include "lapack.h"

#ifdef FORCE_FLOAT
#include "param_float.h"
#include "ldr_float.h"
#include "lm_float.h"
#include "lmm_float.h"//for LMM class,and functions CalcLambda, CalcPve, CalcVgVe
#include "mathfunc_float.h"  //for function CenterVector
#else
#include "param.h"
#include "ldr.h"
#include "lm.h"
#include "bslmm.h"
#include "lmm.h"
#include "mathfunc.h"
#endif

using namespace std;
using namespace Eigen;

void LDR::CopyFromParam (PARAM &cPar)   {
	a_mode=cPar.a_mode;
	d_pace=cPar.d_pace;

	file_bfile=cPar.file_bfile;
	file_geno=cPar.file_geno;
	file_out=cPar.file_out;
	path_out=cPar.path_out;

	ni_total=cPar.ni_total;
	ns_total=cPar.ns_total;
	ni_test=cPar.ni_test;
	ns_test=cPar.ns_test;
	n_cvt=cPar.n_cvt;

// zeng ping LDR
	n_k=cPar.n_k;
	sp=cPar.sp;
	tp=cPar.tp;
	w_step=cPar.w_step;
	s_step=cPar.s_step;
	mixture_no=cPar.mixture_no;
	pD1 =cPar.pD1;
	pD2 =cPar.pD2;
	DIC1=cPar.DIC1;
	BIC1=cPar.BIC1;
	DIC2=cPar.DIC2;
	BIC2=cPar.BIC2;
	LBO=cPar.LBO;
	mnk=cPar.mnk;
	// zeng ping LDR

	indicator_idv=cPar.indicator_idv;
	indicator_snp=cPar.indicator_snp;
	snpInfo=cPar.snpInfo;
	
	return;
}

void LDR::CopyToParam (PARAM &cPar)   {
	/*cPar.time_UtZ=time_UtZ;
	cPar.time_Omega=time_Omega;
	cPar.time_Proposal=time_Proposal;
	cPar.cHyp_initial=cHyp_initial;
	cPar.n_accept=n_accept;
	cPar.pheno_mean=pheno_mean;
	cPar.randseed=randseed;*/
	cPar.LBO=LBO;
	cPar.mnk=mnk;
	cPar.pheno_mean=pheno_mean;
	cPar.mixture_no=mixture_no;

	cPar.pD1 =pD1;
	cPar.pD2 =pD2;
	cPar.DIC1=DIC1;
	cPar.BIC1=BIC1;
	cPar.DIC2=DIC2;
	cPar.BIC2=BIC2;
	return;
}

void LDR::WritevbBeta (const VectorXd Ebeta)  {
	string file_str;
	file_str=path_out+"/"+file_out;
	file_str+=".param.txt";

	ofstream outfile (file_str.c_str(), ofstream::out);
	if (!outfile) {cout<<"error writing file: "<<file_str.c_str()<<endl; return;}

	outfile<<"chr"<<"\t"<<"rs"<<"\t"
	<<"ps"<<"\t"<<"n_miss"<<"\t"<<"b"<<"\t"
	<<"beta"<<"\t"<<"gamma"<<endl;

	size_t t=0;
	for (size_t i=0; i<ns_total; ++i) {
		if (indicator_snp[i]==0) {continue;}	

	outfile<<snpInfo[i].chr<<"\t"<<snpInfo[i].rs_number<<"\t"
	<<snpInfo[i].base_position<<"\t"<<snpInfo[i].n_miss<<"\t";

	outfile<<scientific<<setprecision(6)<<0.0<<"\t";
	outfile<<scientific<<setprecision(6)<<Ebeta(t)<<"\t";
	outfile<<1<<endl;
	t++;
}

	outfile.clear();
	outfile.close();
return;
}


void LDR::WriteParam (const VectorXd alpha) 
{
	string file_str;
	file_str=path_out+"/"+file_out;
	file_str+=".param.txt";

	ofstream outfile (file_str.c_str(), ofstream::out);
	if (!outfile) {cout<<"error writing file: "<<file_str.c_str()<<endl; return;}
	
	outfile<<"chr"<<"\t"<<"rs"<<"\t"
			<<"ps"<<"\t"<<"n_miss"<<"\t"<<"b"<<"\t"
			<<"beta"<<"\t"<<"gamma"<<endl;
	
	size_t t=0;
	for (size_t i=0; i<ns_total; ++i) {
		if (indicator_snp[i]==0) {continue;}		

		outfile<<snpInfo[i].chr<<"\t"<<snpInfo[i].rs_number<<"\t"
				<<snpInfo[i].base_position<<"\t"<<snpInfo[i].n_miss<<"\t";				
		//outfile<<scientific<<setprecision(6)<<gsl_vector_get(alpha, t)<<"\t";
		outfile<<scientific<<setprecision(6)<<alpha(t)<<"\t";
		outfile<<0.0<<"\t"<<0.0<<endl;
		t++;
	}		
	
	outfile.clear();	
	outfile.close();	
	return;
}


void LDR::WriteCoeff (const VectorXd Alpha, const VectorXd Beta) {
	string file_str;
	file_str=path_out+"/"+file_out;
	file_str+=".param.txt";
	
	ofstream outfile (file_str.c_str(), ofstream::out);
	if (!outfile) {cout<<"error writing file: "<<file_str.c_str()<<endl; return;}
	
	outfile<<"chr"<<"\t"<<"rs"<<"\t"<<"ps"<<"\t"<<"n_miss"
	<<"\t"<<"b"<<"\t"<<"beta"<<"\t"<<"gamma"<<endl;
	
	size_t t=0;
	for (size_t i=0; i<ns_total; ++i) {
		if (indicator_snp[i]==0) {continue;}	
		outfile<<snpInfo[i].chr<<"\t"<<snpInfo[i].rs_number<<"\t"
		<<snpInfo[i].base_position<<"\t"<<snpInfo[i].n_miss<<"\t";
		outfile<<scientific<<setprecision(6)<<Alpha(t)<<"\t";
		outfile<<scientific<<setprecision(6)<<Beta(t)<<"\t";
		outfile<<(size_t)((Beta(t)!=0)*1)<<endl;
		t++;
	}

	outfile.clear();
	outfile.close();
	return;
}


void LDR::WritelmmBeta (const VectorXd alpha)  {
	string file_str;
	file_str=path_out+"/"+file_out;
	file_str+=".lmm.param.txt";

	ofstream outfile (file_str.c_str(), ofstream::out);
	if (!outfile) {cout<<"error writing file: "<<file_str.c_str()<<endl; return;}
	
	outfile<<"chr"<<"\t"<<"rs"<<"\t"
			<<"ps"<<"\t"<<"n_miss"<<"\t"<<"b"<<"\t"
			<<"beta"<<"\t"<<"gamma"<<endl;
	
	size_t t=0;
	for (size_t i=0; i<ns_total; ++i) {
		if (indicator_snp[i]==0) {continue;}		

		outfile<<snpInfo[i].chr<<"\t"<<snpInfo[i].rs_number<<"\t"
				<<snpInfo[i].base_position<<"\t"<<snpInfo[i].n_miss<<"\t";				
		//outfile<<scientific<<setprecision(6)<<gsl_vector_get(alpha, t)<<"\t";
		outfile<<scientific<<setprecision(6)<<alpha(t)<<"\t";
		outfile<<0.0<<"\t"<<0.0<<endl;
		t++;
	}		
	
	outfile.clear();	
	outfile.close();	
	return;
}


void LDR::WritevbAlpha (const VectorXd Ealpha,size_t n_j) {
	string file_str;
	file_str=path_out+"/"+file_out;
	file_str+=".alpha.txt";
	
	ofstream outfile (file_str.c_str(), ofstream::out);
	if (!outfile) {cout<<"error writing file: "<<file_str.c_str()<<endl; return;}
	
	for (size_t i=0; i<n_j; i++) {
	outfile<<scientific<<setprecision(6)<<Ealpha(i)<<endl;
	}
	
	outfile.clear();
	outfile.close();
	return;
}

void LDR::WriteLike (const VectorXd Like,size_t t) {
	string file_str;
	file_str=path_out+"/"+file_out;
	file_str+=".like.txt";
	
	ofstream outfile (file_str.c_str(), ofstream::out);
	if (!outfile) {cout<<"error writing file: "<<file_str.c_str()<<endl; return;}
	
	for (size_t i=0; i<t; i++) {
	outfile<<scientific<<setprecision(6)<<Like(i)<<endl;
	}
	
	outfile.clear();
	outfile.close();
	return;
}



void LDR::WritevbVariance (const double variance) {
	string file_str;
	file_str=path_out+"/"+file_out;
	file_str+=".variance.txt";
	
	ofstream outfile (file_str.c_str(), ofstream::out);
	if (!outfile) {cout<<"error writing file: "<<file_str.c_str()<<endl; return;}
	
  outfile<<scientific<<setprecision(6)<<variance<<endl;
  
	outfile.clear();
	outfile.close();
	return;
}

void LDR::WritevbELBO (const VectorXd ELBO,size_t n_int) {
	string file_str;
	file_str=path_out+"/"+file_out;
	file_str+=".ELBO.txt";
	
	ofstream outfile (file_str.c_str(), ofstream::out);
	if (!outfile) {cout<<"error writing file: "<<file_str.c_str()<<endl; return;}
	
		for (size_t i=1; i<=n_int; i++) {
		outfile<<scientific<<setprecision(6)<<i<<"\t"<< ELBO(i)<<endl;
	}
	
	outfile.clear();
	outfile.close();
	return;
}

void LDR::WritegibbsBeta (const VectorXd Ebeta) {
	string file_str;
	file_str=path_out+"/"+file_out;
	file_str+=".param.txt";
	
	ofstream outfile (file_str.c_str(), ofstream::out);
	if (!outfile) {cout<<"error writing file: "<<file_str.c_str()<<endl; return;}
	
	outfile<<"chr"<<"\t"<<"rs"<<"\t"
	<<"ps"<<"\t"<<"n_miss"<<"\t"<<"b"<<"\t"
	<<"beta"<<"\t"<<"gamma"<<endl;
	
	size_t t=0;
	for (size_t i=0; i<ns_total; ++i) {
		if (indicator_snp[i]==0) {continue;}
		
		outfile<<snpInfo[i].chr<<"\t"<<snpInfo[i].rs_number<<"\t"
		<<snpInfo[i].base_position<<"\t"<<snpInfo[i].n_miss<<"\t";
		
		outfile<<scientific<<setprecision(6)<<0.0<<"\t";
		outfile<<scientific<<setprecision(6)<<Ebeta(t)<<"\t";
		outfile<<1<<endl;
		t++;
	}
	
	outfile.clear();
	outfile.close();
	return;
}


void LDR::WritegibbsAlpha (const VectorXd Ealpha, size_t n_j) {
	string file_str;
	file_str=path_out+"/"+file_out;
	file_str+=".alpha.txt";
	
	ofstream outfile (file_str.c_str(), ofstream::out);
	if (!outfile) {cout<<"error writing file: "<<file_str.c_str()<<endl; return;}
	
	for (size_t i=0; i<n_j; i++) {
		outfile<<scientific<<setprecision(6)<<Ealpha(i)<<endl;
	}
	
	outfile.clear();
	outfile.close();
	return;
}


void LDR::WritegibbsVariance (const double variance) {
	string file_str;
	file_str=path_out+"/"+file_out;
	file_str+=".variance.txt";
	
	ofstream outfile (file_str.c_str(), ofstream::out);
	if (!outfile) {cout<<"error writing file: "<<file_str.c_str()<<endl; return;}
	
	outfile<<scientific<<setprecision(6)<<variance<<endl;
	
	outfile.clear();
	outfile.close();
	return;
}


///////////////////////////////////////////
/////// some functions for LDR ///////////
///////////////////////////////////////////
// function 1: E(log(1-vk))
double sum_Elogvl(const VectorXd lambda_k, const VectorXd kappa_k, size_t k) {
	double sumElogvl = 0;
	if (k==0) sumElogvl = 0;
	else {
		for (size_t j = 0; j<k; j++)	{
		sumElogvl += gsl_sf_psi(lambda_k(j))-gsl_sf_psi(kappa_k(j)+lambda_k(j));
		}
	}
	return sumElogvl;
}

// function 2: the first term of lambda_k
double sum_lambda_k(const MatrixXd pik_beta, size_t k, size_t n_k) {
	double slambda_k =0;
	if ((k+1)>(n_k-1)) {slambda_k =0;}
	else {
		for (size_t j=k+1; j<n_k-1; j++) {
		slambda_k += pik_beta.col(j).sum();
		}
	}
	return slambda_k;
}

// function 3: the second term of b_lambda(i.e, sum(Elog(1-vk)))
double sum_b_lambda(const VectorXd lambda_k,const VectorXd kappa_k,size_t n_k) {
	double  sb_lambda=0;
	for (size_t k=0; k<n_k-1; k++) {
		sb_lambda += gsl_sf_psi(lambda_k(k))-gsl_sf_psi(kappa_k(k)+lambda_k(k));
	}
	return sb_lambda;
}

// function 4
double ELBO1(const VectorXd a_k,const VectorXd b_k,size_t n_k) {
	double  sum_ELBO1=0;
	for (size_t k=1; k<n_k; k++) {
		sum_ELBO1 += lgamma(a_k(k))-a_k(k)*log(b_k(k))+a_k(k);
	}
	return sum_ELBO1;
}

// function 5
double ELBO2(const VectorXd kappa_k,const VectorXd lambda_k,size_t n_k) {
	double  sum_ELBO2=0;
	for (size_t k=0; k<n_k-1; k++) {
		sum_ELBO2 += lgamma(kappa_k(k))+lgamma(lambda_k(k))-
		lgamma(kappa_k(k)+lambda_k(k));
	}
	return sum_ELBO2;
}

// function 6
double ELBO3(const MatrixXd pik_beta,const MatrixXd sik2_beta) {
	double  sum_ELBO3=0;
	MatrixXd pik_sik2;
	MatrixXd sik2_betax;
	sik2_betax = sik2_beta*2*3.141593*2.718282;
	pik_sik2=((pik_beta.array()+1e-10).array().log()-0.5*(sik2_betax.array()+1e-10).log());
	sum_ELBO3=(pik_sik2.array()*pik_beta.array()).sum();
	return sum_ELBO3;
}

double ELBO3(const VectorXd pik_beta,const VectorXd sik2_beta) {
	double  sum_ELBO3=0;
	VectorXd pik_sik2;
	VectorXd sik2_betax;
	sik2_betax = sik2_beta*2*3.141593*2.718282;
	pik_sik2=((pik_beta.array()+1e-10).array().log()-0.5*(sik2_betax.array()+1e-10).log());
	sum_ELBO3=(pik_sik2.array()*pik_beta.array()).sum();
	return sum_ELBO3;
}

double ELBO4(const VectorXd pik_beta) {
	double  sum_ELBO4=0;
	VectorXd pik_sik2;
	VectorXd pik_sik3;
	pik_sik2=(pik_beta.array()+1e-10).array().log();
	pik_sik2=pik_beta.array()*pik_sik2.array();
	pik_sik3=(1-pik_beta.array()).array().log();
	pik_sik3=(1-pik_beta.array()).array()*pik_sik3.array();
	sum_ELBO4=(pik_sik2+pik_sik3).sum();
	return sum_ELBO4;
	pik_sik2.resize(0);
	pik_sik3.resize(0);
}

double ELBO5(const VectorXd pik_beta, const VectorXd sik2_beta) {
	double  sum_ELBO5=0;
	VectorXd X;
	X=(sik2_beta*2*3.14).array().log()+1;
	sum_ELBO5=pik_beta.dot(X)/2;
	return sum_ELBO5;
	X.resize(0);
}

//////////////////////////////////////////////
  //////////////////////////////////////////////
double sum_Elogvl2(const VectorXd vk, size_t k) {
	double sumElogvl = 0;
	if (k==0) sumElogvl = 0;
	else {
		for (size_t j = 0; j<k; j++)	{
		sumElogvl += log(1-vk(j));
		}
	}
	return sumElogvl;
}

double log_sigma2b(const VectorXd D, const VectorXd y_res, double sigma2b, double sigma2e, double ae, double be) {
	VectorXd H0 = (sigma2b*D).array()+1;
	double  a = H0.array().log().sum();
		   H0 = y_res.array()/H0.array();
	double  b = H0.dot(y_res);
	double  c = (ae+1)*log(sigma2b)+be/sigma2b;
	double  log_density=-0.5*a-0.5*b/sigma2e-c;
	return  log_density;
}

// log(p(y|theta))
double logLike(const VectorXd D, const VectorXd y_res, double sigma2b, double sigma2e) {
	VectorXd H0 = (sigma2b*D).array()+1;
	double  a = H0.array().log().sum();
		   H0 = y_res.array()/H0.array();
	double  b = H0.dot(y_res);
	size_t  n_idv = y_res.size();
	double  c = n_idv*(log(sigma2e)+log(2*3.1415926));
	double  log_density=-0.5*a-0.5*b/sigma2e-0.5*c;
	return  log_density;
}


double log_h(const VectorXd D, const VectorXd y_res, double h, double sigma2e, double ae, double be) {
	double sigma2b = h/(1-h);
	VectorXd H0 = (sigma2b*D).array()+1;
	double  a = H0.array().log().sum();
		   H0 = y_res.array()/H0.array();
	double  b = H0.dot(y_res);
	double  c = (ae+1)*log(sigma2b)+be/sigma2b;
	double  log_density=-0.5*a-0.5*b/sigma2e-c-2*log(1-h);
	return  log_density;
}


void LDR::call_sig_allsnp_no(gsl_matrix *W, vector<int> *sig_subsnp_no, PARAM &cPar)
{
	//PARAM cPar;
	////////////////////////////Fit LM 
	vector<double> lm_pvalue;
	vector<double> lm_pvalue0;
	//vector<double> lm_beta;
	//vector<double> lm_se;
	//vector<string> lm_rs;

	gsl_matrix *Y=gsl_matrix_alloc (cPar.ni_test, cPar.n_ph);
	cPar.CopyCvtPhen (W, Y, 0);
	if (cPar.n_ph==1) {
	LM cLm;
	cLm.CopyFromParam(cPar);
	gsl_vector_view Y_col=gsl_matrix_column (Y, 0);
	if (!cPar.file_gene.empty()) {
	cLm.AnalyzeGene (W, &Y_col.vector); 
			}
		else if (!cPar.file_bfile.empty()) {
				cLm.AnalyzePlink (W, &Y_col.vector);
				} 
			else if (!cPar.file_oxford.empty()) {
				cLm.Analyzebgen (W, &Y_col.vector);
					} 
				else {
					cLm.AnalyzeBimbam (W, &Y_col.vector);
					}
	cLm.Call_lm_Pvalue(&lm_pvalue);
	//cLm.Call_lm_Beta(&lm_beta);
	//cLm.Call_lm_Se(&lm_se);
	//cLm.Call_lm_Rs(&lm_rs);
	}
		///Significant SNPs screening
	cout<<endl;
	cout<<"Significant SNPs screening..."<<endl;
	lm_pvalue0=lm_pvalue;
	sort(lm_pvalue0.begin(),lm_pvalue0.end());

	double cutoff=cPar.pcutoff;//lm_pvalue.size()
	if (cutoff>lm_pvalue.size())
	{
		cout<<endl;
		cout<<"Due to the option of -m ";cout<<(size_t)cutoff;
		cout<<";";cout<<" selected number of the top marginally "<<endl;
		cout<<"significant SNPs is larger than the number of analyzed SNPs; please "<<endl;
		cout<<"check that, otherwise, all SNPs will be included into the naive MCMC "<<endl;
		cout<<"sampling, which will lead to the increase of the computaional time."<<endl;
		cout<<endl;
		cutoff=2;
	}
	else cutoff=lm_pvalue0[cutoff];//lm_pvalue.size()

	VectorXd  sig_index(lm_pvalue.size());
	for (size_t i=0; i<lm_pvalue.size(); ++i) {
		sig_index(i)=(lm_pvalue[i]<cutoff)*1;
		}
	VectorXd sig_allsnp_no(lm_pvalue.size());
	VectorXd sig_subsnp_no_index(sig_index.sum());
	VectorXd nonsig_subsnp_no(lm_pvalue.size()-sig_index.sum());
	//VectorXd sig_subsnp_se(sig_index.sum());
	VectorXd nonsig_subsnp_se(lm_pvalue.size()-sig_index.sum());
	int icount1=0;
	int icount2=0;
	for (size_t i=0; i<lm_pvalue.size(); ++i) {
		sig_index(i)=(lm_pvalue[i]<cutoff)*1;
		sig_allsnp_no(i)=sig_index(i)*i;
		//cout<<lm_rs[i]<<"\t";
		//cout<<setprecision(6)<<lm_beta[i]<<"\t"<<lm_se[i];
		//cout<<"\t"<<lm_pvalue[i];
		//cout<<setprecision(0)<<"\t"<<sig_index(i);
		//cout<<"\t"<<sig_allsnp_no(i)<<endl;
		if(((lm_pvalue[i]<cutoff)*1)>0) {
			sig_subsnp_no_index(icount1)=sig_allsnp_no(i);
			//sig_subsnp_se(icount1)=lm_se[i];
			icount1++;}
		if(((lm_pvalue[i]<cutoff)*1)==0) {
			nonsig_subsnp_no(icount2)=i;
			//nonsig_subsnp_se(icount2)=lm_se[i];
			icount2++;}
		}
		//cout<<"number of analyzed polygenic SNPs = "<<setprecision(0)<<lm_pvalue.size()-sig_index.sum()<<endl;
		//cout<<"number of analyzed mixture SNPs   = "<<setprecision(0)<<sig_index.sum()<<endl;
		//cout<<setprecision(6)<<endl;
		///Significant SNPs screening
		for (size_t i=0; i<(size_t)sig_subsnp_no_index.size(); ++i) {
			//cout <<sig_subsnp_no_index(i)<<endl;
			//(*sig_subsnp_no)(i) = sig_subsnp_no_index(i);
			(*sig_subsnp_no).push_back(sig_subsnp_no_index(i));
			}
	sig_index.resize(0);
	sig_allsnp_no.resize(0);
	sig_subsnp_no_index.resize(0);
	nonsig_subsnp_no.resize(0);
	//sig_subsnp_se.resize(0);
	nonsig_subsnp_se.resize(0);
	gsl_matrix_free(Y);
}


void LDR::call_sig_index(gsl_matrix *W, vector<int> *sig_index, PARAM &cPar)
{
	//PARAM cPar;
	////////////////////////////Fit LM 
	vector<double> lm_pvalue;
	vector<double> lm_beta;
	vector<double> lm_se;
	vector<string> lm_rs;
	gsl_matrix *Y=gsl_matrix_alloc (cPar.ni_test, cPar.n_ph);
	cPar.CopyCvtPhen (W, Y, 0);
	if (cPar.n_ph==1) {
	LM cLm;
	cLm.CopyFromParam(cPar);
	gsl_vector_view Y_col=gsl_matrix_column (Y, 0);
	if (!cPar.file_gene.empty()) {
	cLm.AnalyzeGene (W, &Y_col.vector); 
			}
		else if (!cPar.file_bfile.empty()) {
				cLm.AnalyzePlink (W, &Y_col.vector);
				} 
			else if (!cPar.file_oxford.empty()) {
				cLm.Analyzebgen (W, &Y_col.vector);
					} 
				else {
					cLm.AnalyzeBimbam (W, &Y_col.vector);
					}
	cLm.Call_lm_Pvalue(&lm_pvalue);
	cLm.Call_lm_Beta(&lm_beta);
	cLm.Call_lm_Se(&lm_se);
	cLm.Call_lm_Rs(&lm_rs);
	}
		///Significant SNPs screening
	double cutoff=cPar.pcutoff;//lm_pvalue.size()
	VectorXd  sig_index0(lm_pvalue.size());
	for (size_t i=0; i<lm_pvalue.size(); ++i) {
		sig_index0(i)=(lm_pvalue[i]<cutoff)*1;
		}

	VectorXd sig_allsnp_no(lm_pvalue.size());
	VectorXd sig_subsnp_no_index(sig_index0.sum());
	VectorXd nonsig_subsnp_no(lm_pvalue.size()-sig_index0.sum());
	VectorXd sig_subsnp_se(sig_index0.sum());
	VectorXd nonsig_subsnp_se(lm_pvalue.size()-sig_index0.sum());
	int icount1=0;
	int icount2=0;
	for (size_t i=0; i<lm_pvalue.size(); ++i) {
		sig_index0(i)=(lm_pvalue[i]<cutoff)*1;
		sig_allsnp_no(i)=sig_index0(i)*i;
		//cout<<lm_rs[i]<<"\t";
		//cout<<setprecision(6)<<lm_beta[i]<<"\t"<<lm_se[i];
		//cout<<"\t"<<lm_pvalue[i];
		//cout<<setprecision(0)<<"\t"<<sig_index0(i);
		//cout<<"\t"<<sig_allsnp_no(i)<<endl;
		if(((lm_pvalue[i]<cutoff)*1)>0) {
			sig_subsnp_no_index(icount1)=sig_allsnp_no(i);
			sig_subsnp_se(icount1)=lm_se[i];
			icount1++;}
		if(((lm_pvalue[i]<cutoff)*1)==0) {
			nonsig_subsnp_no(icount2)=i;
			nonsig_subsnp_se(icount2)=lm_se[i];
			icount2++;}
		}
	///Significant SNPs screening
	for (size_t i=0; i<(size_t)sig_index0.size(); ++i) {
		//(*sig_index)(i) = sig_index0(i);
		(*sig_index).push_back(sig_index0(i));
		}

	sig_index0.resize(0);
	sig_allsnp_no.resize(0);
	sig_subsnp_no_index.resize(0);
	nonsig_subsnp_no.resize(0);
	sig_subsnp_se.resize(0);
	nonsig_subsnp_se.resize(0);
	gsl_matrix_free(Y);
}

void LDR::call_Chr(gsl_matrix *W, vector<int> *Chr, PARAM &cPar)
{
	//PARAM cPar;
	////////////////////////////Fit LM 
	//vector<double> lm_pvalue;
	vector<int> lm_chr;
	//vector<double> lm_beta;
	//vector<double> lm_se;
	//vector<string> lm_rs;

	gsl_matrix *Y=gsl_matrix_alloc (cPar.ni_test, cPar.n_ph);
	cPar.CopyCvtPhen (W, Y, 0);
	if (cPar.n_ph==1) {
	LM cLm;
	cLm.CopyFromParam(cPar);
	gsl_vector_view Y_col=gsl_matrix_column (Y, 0);
	if (!cPar.file_gene.empty()) {
	cLm.AnalyzeGene (W, &Y_col.vector); 
			}
		else if (!cPar.file_bfile.empty()) {
				cLm.AnalyzePlink (W, &Y_col.vector);
				} 
			else if (!cPar.file_oxford.empty()) {
				cLm.Analyzebgen (W, &Y_col.vector);
					} 
				else {
					cLm.AnalyzeBimbam (W, &Y_col.vector);
					}
	//cLm.Call_lm_Pvalue(&lm_pvalue);
	//cLm.Call_lm_Beta(&lm_beta);
	//cLm.Call_lm_Se(&lm_se);
	//cLm.Call_lm_Rs(&lm_rs);
	cLm.Call_lm_Chr(&lm_chr);
	}

	for (size_t i=0; i<(size_t)lm_chr.size(); ++i) {
		(*Chr).push_back(lm_chr[i]);
		}
	gsl_matrix_free(Y);
}

void LDR::call_Position(gsl_matrix *W, vector<int> *ps, PARAM &cPar)
{
	//PARAM cPar;
	////////////////////////////Fit LM 
	//vector<double> lm_pvalue;
	vector<int> lm_ps;
	//vector<double> lm_beta;
	//vector<double> lm_se;
	//vector<string> lm_rs;

	gsl_matrix *Y=gsl_matrix_alloc (cPar.ni_test, cPar.n_ph);
	cPar.CopyCvtPhen (W, Y, 0);
	if (cPar.n_ph==1) {
	LM cLm;
	cLm.CopyFromParam(cPar);
	gsl_vector_view Y_col=gsl_matrix_column (Y, 0);
	if (!cPar.file_gene.empty()) {
	cLm.AnalyzeGene (W, &Y_col.vector); 
			}
		else if (!cPar.file_bfile.empty()) {
				cLm.AnalyzePlink (W, &Y_col.vector);
				} 
			else if (!cPar.file_oxford.empty()) {
				cLm.Analyzebgen (W, &Y_col.vector);
					} 
				else {
					cLm.AnalyzeBimbam (W, &Y_col.vector);
					}
	//cLm.Call_lm_Pvalue(&lm_pvalue);
	//cLm.Call_lm_Beta(&lm_beta);
	//cLm.Call_lm_Se(&lm_se);
	//cLm.Call_lm_Rs(&lm_rs);
	cLm.Call_lm_Position(&lm_ps);
	}

	for (size_t i=0; i<(size_t)lm_ps.size(); ++i) {
		(*ps).push_back(lm_ps[i]);
		}
	gsl_matrix_free(Y);
}


void LDR::call_Chr_Position(gsl_matrix *W, vector<int> *Chr, vector<int> *ps, PARAM &cPar)
{
	//PARAM cPar;
	////////////////////////////Fit LM 
	//vector<double> lm_pvalue;
	vector<int> lm_chr;
	vector<int> lm_ps;
	//vector<double> lm_beta;
	//vector<double> lm_se;
	//vector<string> lm_rs;

	gsl_matrix *Y=gsl_matrix_alloc (cPar.ni_test, cPar.n_ph);
	cPar.CopyCvtPhen (W, Y, 0);
	if (cPar.n_ph==1) {
	LM cLm;
	cLm.CopyFromParam(cPar);
	gsl_vector_view Y_col=gsl_matrix_column (Y, 0);
	if (!cPar.file_gene.empty()) {
	cLm.AnalyzeGene (W, &Y_col.vector); 
			}
		else if (!cPar.file_bfile.empty()) {
				cLm.AnalyzePlink (W, &Y_col.vector);
				} 
			else if (!cPar.file_oxford.empty()) {
				cLm.Analyzebgen (W, &Y_col.vector);
					} 
				else {
					cLm.AnalyzeBimbam (W, &Y_col.vector);
					}
	//cLm.Call_lm_Pvalue(&lm_pvalue);
	//cLm.Call_lm_Beta(&lm_beta);
	//cLm.Call_lm_Se(&lm_se);
	//cLm.Call_lm_Rs(&lm_rs);
	cLm.Call_lm_Chr_Position(&lm_chr, &lm_ps);
	//cLm.Call_lm_Position(&lm_ps);
	}

	for (size_t i=0; i<(size_t)lm_chr.size(); ++i) {
		(*Chr).push_back(lm_chr[i]);
		 (*ps).push_back(lm_ps[i]);
		}
	gsl_matrix_free(Y);
}


void LDR::call_Pvalue(gsl_matrix *W, vector<double> *pvalue, PARAM &cPar)
{
	//PARAM cPar;
	////////////////////////////Fit LM 
	vector<double> lm_pvalue;
	//vector<double> lm_ps;
	//vector<double> lm_beta;
	//vector<double> lm_se;
	//vector<string> lm_rs;

	gsl_matrix *Y=gsl_matrix_alloc (cPar.ni_test, cPar.n_ph);
	cPar.CopyCvtPhen (W, Y, 0);
	if (cPar.n_ph==1) {
	LM cLm;
	cLm.CopyFromParam(cPar);
	gsl_vector_view Y_col=gsl_matrix_column (Y, 0);
	if (!cPar.file_gene.empty()) {
	cLm.AnalyzeGene (W, &Y_col.vector); 
			}
		else if (!cPar.file_bfile.empty()) {
				cLm.AnalyzePlink (W, &Y_col.vector);
				} 
			else if (!cPar.file_oxford.empty()) {
				cLm.Analyzebgen (W, &Y_col.vector);
					} 
				else {
					cLm.AnalyzeBimbam (W, &Y_col.vector);
					}
	cLm.Call_lm_Pvalue(&lm_pvalue);
	//cLm.Call_lm_Beta(&lm_beta);
	//cLm.Call_lm_Se(&lm_se);
	//cLm.Call_lm_Rs(&lm_rs);
	//cLm.Call_lm_Position(&lm_ps);
	}

	for (size_t i=0; i<(size_t)lm_pvalue.size(); ++i) {
		(*pvalue).push_back(lm_pvalue[i]);
		}
	gsl_matrix_free(Y);
}


void LDR::LDscreen(gsl_matrix *LDblock0, MatrixXd LDblock, gsl_matrix *W, gsl_vector *y, gsl_matrix *UtX, vector<int> *selectX_snp, PARAM &cPar)
{
	//cout<<endl<<endl;
	cout<<endl<<"1) Obtain LD block matrix"<<endl;
	for (size_t i=0; i<LDblock0->size1; ++i)  {
		for (size_t j=0; j<3; ++j) {
			LDblock(i,j) = gsl_matrix_get(LDblock0,i,j);
		}
		//cout<<i+1<<" "<<LDblock.row(i)<<endl;
	}

	vector<int> XChr, Xps;
	cout<<endl<<"2) Obtain Chromosome and SNP-Position from the present data"<<endl;
	call_Chr_Position (W, &XChr, &Xps, cPar);
	MatrixXi Chr_ps  = MatrixXi::Zero(XChr.size(),2);
	for (size_t i=0; i<XChr.size(); ++i) {
		Chr_ps(i,0) = XChr[i];
		Chr_ps(i,1) = Xps[i];
		//cout<<i+1<<" "<<Chr_ps.row(i)<<endl;
	}

	VectorXi index  = VectorXi::Zero((size_t)Chr_ps.rows());
	size_t   icount = 0;
	for (size_t j=0; j<(size_t)LDblock.rows(); ++j) {
		for (size_t i=0; i<(size_t)Chr_ps.rows(); ++i) {
			if ((size_t)Chr_ps(i,0) == (size_t)LDblock(j,0)) {
				if ((Chr_ps(i,1)>=LDblock(j,1)) && (Chr_ps(i,1)<=LDblock(j,2))) {
					index (icount) = j+1;
					icount++;
				}
				if (icount>=(size_t)Chr_ps.rows()) {
					icount=(size_t)Chr_ps.rows()-1;
					}
			}
		//cout<<j+1<<" "<<i+1<<" "<<icount<<" "<<index(icount)<<endl;
		}
	}

	cout<<endl<<"3) Specify SNPs into LD blocks, note some blocks may have no SNPs"<<endl;
	size_t    total_ld = index.maxCoeff();
	//total_ld is the total number of ld block in the present data
	MatrixXi indexSNP1  = MatrixXi::Zero(index.maxCoeff(), 3);
	for (size_t j=0; j<total_ld; ++j) {
		for (size_t i=0; i<XChr.size(); ++i) {
			if ((size_t)index(i)==(j+1)) {
				indexSNP1(j,0) += 1;
				}
		}
		if (j==0) {
			indexSNP1(j,2) = indexSNP1(j,0);
			indexSNP1(j,1) = 1 ;
			}
		else      {
			indexSNP1(j,2) = indexSNP1(j,0) + indexSNP1(j-1,2);
			indexSNP1(j,1) = indexSNP1(j-1,2) + 1;
			}
			//cout<<j+1<<" "<<indexSNP1.row(j)<<endl;
	}

	cout<<endl<<"4) Romove LD blocks having no SNPs, ";
	vector<int> index0, index1, index2;
	for (size_t i=0; i<(size_t)indexSNP1.rows(); ++i) {
		if (indexSNP1(i,0)>0) {
			index0.push_back(indexSNP1(i,0));
			index1.push_back(indexSNP1(i,1));
			index2.push_back(indexSNP1(i,2));
		}
	}

	MatrixXi indexSNP2(index0.size(),3);
	for (size_t i=0; i<index0.size(); ++i) {
		indexSNP2(i,0)=index0[i];
		indexSNP2(i,1)=index1[i];
		indexSNP2(i,2)=index2[i];
		//cout<<i+1<<" "<<indexSNP2.row(i)<<endl;
	}

	MatrixXi indexSNP;
	if (indexSNP2(indexSNP2.rows()-1,0)>1) {indexSNP=indexSNP2;}
	if (indexSNP2(indexSNP2.rows()-1,0)==1) {
		indexSNP=indexSNP2.topRows(indexSNP2.rows()-1);
		indexSNP(indexSNP.rows()-1,0)+= 1;
		indexSNP(indexSNP.rows()-1,2) = indexSNP2(indexSNP2.rows()-1,2);
	}

	for (size_t i=0; i<(size_t)indexSNP.rows(); ++i) {
		cout<<i+1<<" "<<indexSNP.row(i)<<endl;
	}

	cout<<"after this now"<<" "<<indexSNP.rows()<<" "<<"LD blocks are left"<<endl;
	cout<<endl<<"5) Transform SNP matrix from gsl to eigen"<<" ";
	cout<<"and center phenotype"<<endl;
	size_t n  = y->size, p  = UtX->size2;
	MatrixXd G0  (n, p); VectorXd y00 (n);
	for (size_t i=0; i<n; i++) {
		for (size_t j=0; j<p; j++) {
			G0(i,j) = gsl_matrix_get(UtX,i,j);
		}
		y00(i) = gsl_vector_get(y,i);
		if (i%((size_t)(n*0.2))==0 || i%(n-1)==0) {
			ProgressBar ("   Transform  ",i,n-1);}
	}
	y00 = y00.array()-y00.mean();

	VectorXd x_col(n);
	vector<int> select_block, select_snp;
	//vector<int>    select_snp;
	vector<double> select_pvalue;

	total_ld = indexSNP.rows();
	cout<<endl<<endl;
	cout<<"6) Now conduct conditional step-wise regression screen"<<endl;
	for (size_t j=0; j<total_ld; ++j) {
		size_t L = indexSNP(j,0);
		MatrixXd Gsub    = G0.middleCols(indexSNP(j,1)-1, L);
		MatrixXd p_block = MatrixXd::Ones(L,1);

		VectorXd beta    = VectorXd::Ones(L);
		VectorXd se_wald = VectorXd::Ones(L);

		size_t t  = 0;
		//double df = n-1;
		//VectorXd cutoff (5);
		//cutoff << 0.05, 0.05, 0.05, 0.05, 0.05;
		double df=n-1, cutoff=cPar.pcutoff, min_p=0, sigma2;
		VectorXd yx = y00;
		while ((min_p < cutoff) && (t<5)) {
			ptrdiff_t mi = 0, mj = 0;
			for (size_t m = 0; m < L; ++m) {
				x_col  = Gsub.col(m);
				beta(m)= 1/(x_col.dot(x_col))*(x_col.dot(yx));
				sigma2 = 1/df*(yx-x_col*beta(m)).dot(yx-x_col*beta(m));
				se_wald(m) = 1/(x_col.dot(x_col))*sigma2;
				p_block(m,0)= gsl_cdf_fdist_Q (beta(m)*beta(m)/se_wald(m),1.0,df);
				}

			min_p   = p_block.minCoeff(&mi,&mj);

			if (j==0) {
				//cout<<j+1<<" "<<mi<<" "<<mj;
				//cout<<" "<<setprecision(8)<<min_p<<endl;
				select_block.push_back(j+1);
				select_snp.push_back(mi);
				select_pvalue.push_back(min_p);
			}
			else {
				//cout<<j+1<<" "<<mi+indexSNP(j-1,2)<<" "<<mj;
				//cout<<" "<<setprecision(8)<<min_p<<endl;
				select_block.push_back(j+1);
				select_snp.push_back(mi+indexSNP(j-1,2));
				select_pvalue.push_back(min_p);
			}
			yx      = yx-Gsub.col(mi)*beta(mi);
			t++;
		}

		if (j%((size_t)(total_ld*0.2))==0 || j%(total_ld-1)==0) {ProgressBar ("   Conditional screen  ",j,total_ld-1);}
	}

	G0.resize(0,0);
	cout<<endl<<endl<<"7) Remove SNPs with p-value greater than the given cut-off";
	cout<<" "<<setprecision(4)<<cPar.pcutoff<<endl;
	for (size_t i = 0; i < select_snp.size(); ++i) {
		if (select_pvalue[i]<cPar.pcutoff) {
			//cout<<select_block[i] <<" ";
			//cout<<select_snp[i]   <<" ";
			//cout<<setprecision(8)<<select_pvalue[i]<<" ";
			//cout<<endl;
			(*selectX_snp).push_back(select_snp[i]);
		}
	}

	XChr.clear();
	Xps.clear();
	Chr_ps.resize(0,0);
	index.resize(0);
	indexSNP.resize(0,0);
	y00.resize(0);
	select_block.clear();
	select_snp.clear();
	select_pvalue.clear();
}



void LDR::LDscreen(gsl_matrix *LDblock0, MatrixXd LDblock, gsl_matrix *W, gsl_vector *y, MatrixXd UtX, vector<int> *selectX_snp, PARAM &cPar)
{
	//cout<<endl<<endl;
	cout<<endl<<"1) Obtain LD block matrix"<<endl;
	for (size_t i=0; i<LDblock0->size1; ++i)  {
		for (size_t j=0; j<3; ++j) {
			LDblock(i,j) = gsl_matrix_get(LDblock0,i,j);
		}
		//cout<<i+1<<" "<<LDblock.row(i)<<endl;
	}

	vector<int> XChr, Xps;
	cout<<endl<<"2) Obtain Chromosome and SNP-Position from the present data"<<endl;
	call_Chr_Position (W, &XChr, &Xps, cPar);
	MatrixXi Chr_ps  = MatrixXi::Zero(XChr.size(),2);
	for (size_t i=0; i<XChr.size(); ++i) {
		Chr_ps(i,0) = XChr[i];
		Chr_ps(i,1) = Xps[i];
		//cout<<i+1<<" "<<Chr_ps.row(i)<<endl;
	}

	VectorXi index  = VectorXi::Zero((size_t)Chr_ps.rows());
	size_t   icount = 0;
	for (size_t j=0; j<(size_t)LDblock.rows(); ++j) {
		for (size_t i=0; i<(size_t)Chr_ps.rows(); ++i) {
			if ((size_t)Chr_ps(i,0) == (size_t)LDblock(j,0)) {
				if ((Chr_ps(i,1)>=LDblock(j,1)) && (Chr_ps(i,1)<=LDblock(j,2))) {
					index (icount) = j+1;
					icount++;
				}
				if (icount>=(size_t)Chr_ps.rows()) {
					icount=(size_t)Chr_ps.rows()-1;
					}
			}
		//cout<<j+1<<" "<<i+1<<" "<<icount<<" "<<index(icount)<<endl;
		}
	}

	cout<<endl<<"3) Specify SNPs into LD blocks, note some blocks may have no SNPs"<<endl;
	size_t    total_ld = index.maxCoeff();
	//total_ld is the total number of ld block in the present data
	MatrixXi indexSNP1  = MatrixXi::Zero(index.maxCoeff(), 3);
	for (size_t j=0; j<total_ld; ++j) {
		for (size_t i=0; i<XChr.size(); ++i) {
			if ((size_t)index(i)==(j+1)) {
				indexSNP1(j,0) += 1;
				}
		}
		if (j==0) {
			indexSNP1(j,2) = indexSNP1(j,0);
			indexSNP1(j,1) = 1 ;
			}
		else      {
			indexSNP1(j,2) = indexSNP1(j,0) + indexSNP1(j-1,2);
			indexSNP1(j,1) = indexSNP1(j-1,2) + 1;
			}
			//cout<<j+1<<" "<<indexSNP1.row(j)<<endl;
	}

	cout<<endl<<"4) Romove LD blocks having no SNPs, ";
	vector<int> index0, index1, index2;
	for (size_t i=0; i<(size_t)indexSNP1.rows(); ++i) {
		if (indexSNP1(i,0)>0) {
			index0.push_back(indexSNP1(i,0));
			index1.push_back(indexSNP1(i,1));
			index2.push_back(indexSNP1(i,2));
		}
	}

	MatrixXi indexSNP2(index0.size(),3);
	for (size_t i=0; i<index0.size(); ++i) {
		indexSNP2(i,0)=index0[i];
		indexSNP2(i,1)=index1[i];
		indexSNP2(i,2)=index2[i];
		//cout<<i+1<<" "<<indexSNP2.row(i)<<endl;
	}

	MatrixXi indexSNP;
	if (indexSNP2(indexSNP2.rows()-1,0)>1) {indexSNP=indexSNP2;}
	if (indexSNP2(indexSNP2.rows()-1,0)==1) {
		indexSNP=indexSNP2.topRows(indexSNP2.rows()-1);
		indexSNP(indexSNP.rows()-1,0)+= 1;
		indexSNP(indexSNP.rows()-1,2) = indexSNP2(indexSNP2.rows()-1,2);
	}

	cout<<"after this now"<<" "<<indexSNP.rows()<<" "<<"LD blocks are left"<<endl;
	cout<<endl<<"5) Center phenotype"<<endl;
	size_t n  = y->size;
	VectorXd y00 (n);
	for (size_t i=0; i<n; i++) {
		y00(i) = gsl_vector_get(y,i);
	}
	y00 = y00.array()-y00.mean();

	VectorXd x_col(n);
	vector<int> select_block, select_snp;
	vector<double> select_pvalue;

	total_ld = indexSNP.rows();
	cout<<endl<<"6) Now conduct conditional step-wise regression screen"<<endl;
	for (size_t j=0; j<total_ld; ++j) {
		size_t L = indexSNP(j,0);
		MatrixXd Gsub    = UtX.middleCols(indexSNP(j,1)-1, L);
		MatrixXd p_block = MatrixXd::Ones(L,1);

		VectorXd beta    = VectorXd::Ones(L);
		VectorXd se_wald = VectorXd::Ones(L);

		size_t t  = 0;
		//double df = n-1;
		//VectorXd cutoff (5);
		//cutoff << 0.05, 0.05, 0.05, 0.05, 0.05;
		double df=n-1, cutoff=cPar.pcutoff, min_p=0, sigma2;
		VectorXd yx = y00;
		while ((min_p < cutoff) && (t<5)) {
			ptrdiff_t mi = 0, mj = 0;
			for (size_t m = 0; m < L; ++m) {
				x_col  = Gsub.col(m);
				beta(m)= 1/(x_col.dot(x_col))*(x_col.dot(yx));
				sigma2 = 1/df*(yx-x_col*beta(m)).dot(yx-x_col*beta(m));
				se_wald(m) = 1/(x_col.dot(x_col))*sigma2;
				p_block(m,0)= gsl_cdf_fdist_Q (beta(m)*beta(m)/se_wald(m),1.0,df);
				}

			min_p   = p_block.minCoeff(&mi,&mj);

			if (j==0) {
				//cout<<j+1<<" "<<mi<<" "<<mj;
				//cout<<" "<<setprecision(8)<<min_p<<endl;
				select_block.push_back(j+1);
				select_snp.push_back(mi);
				select_pvalue.push_back(min_p);
			}
			else {
				//cout<<j+1<<" "<<mi+indexSNP(j-1,2)<<" "<<mj;
				//cout<<" "<<setprecision(8)<<min_p<<endl;
				select_block.push_back(j+1);
				select_snp.push_back(mi+indexSNP(j-1,2));
				select_pvalue.push_back(min_p);
			}
			yx      = yx-Gsub.col(mi)*beta(mi);
			t++;
		}

		if (j%((size_t)(total_ld*0.2))==0 || j%(total_ld-1)==0) {ProgressBar ("   Conditional screen  ",j,total_ld-1);}
	}

	//G0.resize(0,0);
	cout<<endl<<endl<<"7) Remove SNPs with p-value greater than the given cut-off";
	cout<<" "<<setprecision(4)<<cPar.pcutoff<<endl;
	for (size_t i = 0; i < select_snp.size(); ++i) {
		if (select_pvalue[i]<cPar.pcutoff) {
			//cout<<select_block[i] <<" ";
			//cout<<select_snp[i]   <<" ";
			//cout<<setprecision(8)<<select_pvalue[i]<<" ";
			//cout<<endl;
			(*selectX_snp).push_back(select_snp[i]);
		}
	}

	XChr.clear();
	Xps.clear();
	Chr_ps.resize(0,0);
	index.resize(0);
	indexSNP.resize(0,0);
	y00.resize(0);
	select_block.clear();
	select_snp.clear();
	select_pvalue.clear();
}

//void LDR::initial_lmm(gsl_matrix *W, gsl_vector *y, gsl_vector *Uty, gsl_matrix *U, gsl_matrix *UtX, gsl_matrix *UtW, gsl_vector *eval, vector<int> *selectX_block, vector<int> *selectX_snp, vector<double> *selectX_pvalue, PARAM &cPar)
//{
//	gsl_matrix *B=gsl_matrix_alloc (y->size, W->size2);
//	gsl_vector *Wbeta=gsl_vector_alloc (W->size2);
//	gsl_vector *se_Wbeta=gsl_vector_alloc (W->size2);
//	gsl_matrix *se_B=gsl_matrix_alloc (y->size, W->size2);
//	gsl_vector *beta  =gsl_vector_alloc (UtX->size2);
//	gsl_vector *H_eval=gsl_vector_alloc (Uty->size);
//	gsl_vector *bv=gsl_vector_alloc (Uty->size);
//	gsl_vector *muw=gsl_vector_alloc (cPar.ni_test);
//	gsl_vector *y_res_w=gsl_vector_alloc (cPar.ni_test);		
//	
//	CalcLambda ('R', eval, UtW, Uty, cPar.l_min, cPar.l_max, cPar.n_region, cPar.l_remle_null, cPar.logl_remle_H0);
//	CalcPve (eval, UtW, Uty, cPar.l_remle_null, cPar.trace_G, cPar.pve_null, cPar.pve_se_null);
//
//	gsl_vector_view xbeta=gsl_matrix_row (B, 0);
//	gsl_vector_view se_xbeta=gsl_matrix_row (se_B, 0);
//	
//	CalcLmmVgVeBeta (eval, UtW, Uty, cPar.l_remle_null, cPar.vg_remle_null, cPar.ve_remle_null, &xbeta.vector, &se_xbeta.vector);
//	
//	for (size_t i=0; i<B->size2; i++) {
//		cPar.beta_remle_null.push_back(gsl_matrix_get(B, 0, i) );
//		gsl_vector_set(Wbeta,i, gsl_matrix_get(B, 0, i));
//		cPar.se_beta_remle_null.push_back(gsl_matrix_get(se_B, 0, i) );
//		gsl_vector_set(se_Wbeta,i, gsl_matrix_get(se_B, 0, i));
//		}
//
//	gsl_vector_memcpy (y_res_w, y);
//	gsl_blas_dgemv (CblasNoTrans, -1, W, Wbeta, 0, muw);
//	gsl_vector_add (y_res_w, muw);
//	CalcUtX (U, y_res_w, Uty);
//	gsl_vector_memcpy (H_eval, eval);
//	gsl_vector_scale (H_eval, cPar.l_remle_null);
//	gsl_vector_add_constant (H_eval, 1.0);
//	gsl_vector_memcpy (bv, Uty);
//	gsl_vector_div (bv, H_eval);	
//	gsl_blas_dgemv (CblasTrans, cPar.l_remle_null/(double)UtX->size2, UtX, bv, 0.0, beta);
//
//	double lambda = cPar.l_remle_null;
//	cout<< endl;
//	cout<<"Estimate for the standard BLUP model" <<endl;
//	cPar.PrintSummary();
//	cout<<"Lambda estimate = " << lambda <<endl;
//}
//
//
//


void Eigen_matrix_get_row (const gsl_matrix *X, const size_t i_col, VectorXd &x_col)
{
  if (i_col < X->size2) {
    for (size_t j=0; j < (size_t)x_col.size(); j++) {
      x_col(j)=gsl_matrix_get(X,j,i_col);
    }
  } else {
    std::cerr << "Error return genotype vector...\n";
    exit(1);
  }
}

//permutation the snp label
VectorXd permute(size_t n_snp)  {
	VectorXd snp_label (n_snp);
	long int randseed;
	gsl_rng *gsl_r;
	time_t rawtime;
	time (&rawtime);
	tm * ptm = gmtime (&rawtime);
	randseed = (unsigned) (ptm->tm_hour%24*3600+ptm->tm_min*60+ptm->tm_sec);
	const gsl_rng_type * gslType;
	gslType = gsl_rng_default; 
	gsl_r = gsl_rng_alloc(gslType); 
	gsl_rng_set(gsl_r, randseed);
	gsl_permutation * perm = gsl_permutation_alloc (n_snp);
	gsl_permutation_init (perm);
	gsl_ran_shuffle (gsl_r, perm->data, n_snp, sizeof(size_t));
	for (size_t i = 0; i<n_snp; i++)	{
		snp_label(i) = perm->data[i];
	}
	return snp_label;
}


void LDR::VB (const MatrixXd UtX, const VectorXd Uty, const MatrixXd UtW,const VectorXd D, const VectorXd Wbeta, const VectorXd se_Wbeta, const VectorXd beta, double lambda) 
{
	
	setprecision(5); 
	clock_t time_begin  = clock();
	size_t n_snp = UtX.cols();
	size_t n_idv = Uty.size();
	size_t n_j   = UtW.cols();
	VectorXd x_col(n_idv);

	VectorXd xtx(n_snp); 
	VectorXd xty(n_snp); 
	VectorXd wtw(n_j); 
	VectorXd wty(n_j);
	VectorXd Utu = VectorXd::Zero(n_idv);
	VectorXd Ue(n_idv);
	VectorXd Ub(n_idv);
	VectorXd bv = VectorXd::Zero(n_idv); //for computation of alpha
	VectorXd V = VectorXd::Zero (n_idv);//compute Utu 
	VectorXd M = VectorXd::Zero (n_idv);//compute Utu
	VectorXd snp_label(n_snp);

	VectorXd Ealpha      = Wbeta;//intercept
	VectorXd post_Ealpha = VectorXd::Zero(n_j);//save intercept
	VectorXd m_alpha     = Wbeta;
	VectorXd s2_alpha    = se_Wbeta.array()*se_Wbeta.array();

	VectorXd WEalpha(n_idv); 
	VectorXd XEbeta(n_idv); //G*beta
	VectorXd Ebeta(n_snp);
	VectorXd post_Ebeta(n_snp);//save beta
	MatrixXd Ebeta2k;
	MatrixXd B1;
	MatrixXd mik_beta = MatrixXd::Zero(n_snp, n_k);

	for (size_t i=0;i<n_snp; i++)  {
		//Eigen_matrix_get_row (UtX, i, x_col);
		x_col  = UtX.col(i);
		xtx(i) = x_col.dot(x_col);
		xty(i) = x_col.dot(Uty);
	}

	for (size_t j=0;j<n_j; j++)  {
		wtw(j) = UtW.col(j).dot(UtW.col(j));
		wty(j) = UtW.col(j).dot(Uty);
	}

	//// initial values for SNP and each normal
	//// component has the same initial values
	for (size_t i=0;i<n_snp; i++)  {
		for (size_t k=1;k<n_k; k++)  {
			mik_beta(i,k) = beta(i);
	}
	}

	mik_beta.col(0).fill(0);
	MatrixXd sik2_beta  =(MatrixXd::Random(n_snp, n_k)).array().abs();
	MatrixXd pik_beta   = MatrixXd::Zero  (n_snp, n_k);
	pik_beta = pik_beta.array() + (double)1/n_k;

	//MatrixXd beta_beta  = mik_beta;
	MatrixXd beta_beta  = MatrixXd::Zero(n_snp, n_k);
	MatrixXd gamma_beta = MatrixXd::Zero(n_snp, n_k);
	MatrixXd post_gamma_beta = gamma_beta;
	MatrixXd post_pik_beta   = gamma_beta;
	gamma_beta = gamma_beta.array() + (double)1/n_k;

	VectorXd sigma2k (n_k);
	VectorXd vk (n_k);
	VectorXd Elogsigmak (n_k);
	VectorXd Elogvk (n_k);
	VectorXd sumElogvl (n_k);

	VectorXd pikexp1(n_k);
	VectorXd pikexp2(n_k);
	VectorXd index;
	VectorXd a_k = (VectorXd::Random(n_k)).array().abs();
	VectorXd b_k = (VectorXd::Random(n_k)).array().abs();
	VectorXd  kappa_k = pik_beta.colwise().sum();
	VectorXd lambda_k = (VectorXd::Random(n_k)).array().abs();

	XEbeta.setZero();//set to zero
	for (size_t i=0; i<n_snp; i++)   {
		x_col   =  UtX.col(i);
		XEbeta += x_col*beta(i);
		}

	WEalpha.setZero(); //set to zero
	for (size_t j=0;j<n_j; j++)   {
		WEalpha += UtW.col(j)*m_alpha(j);
		}

	VectorXd y_res  = Uty - WEalpha - XEbeta;
	double sigma2e0 = ((y_res).dot(y_res))/(n_idv-1);
	//cout << endl << "sigma2e0 = " << sigma2e0 << endl;
	cout << endl;
	double a_e  = (2*n_idv + n_snp)/2;
	double b_e  = (a_e-1)*sigma2e0;

	double ak   = 21;
	if (lambda < 0.01)     {lambda =    0.01;}
	else if (lambda > 100) {lambda =     100;}
	else                   {lambda =  lambda;}
	double bk0  = lambda*(ak-1)/n_snp;

	VectorXd bk(n_k);
	bk(0) = bk0;
	for (size_t i=1;i<n_k; i++)  {
		//bk(i) = bk(i-1)*3.7*sqrt(pow(i,i));
		bk(i) = bk(i-1)*1.7*sqrt(pow(i,i));
		//cout<<bk(i)/(ak-1)*sigma2e0<<" ";
		}

	for (size_t i=0;i<n_snp; i++)  {
		for (size_t k=0;k<n_k; k++)  {
			sik2_beta(i,k) = bk(k)/(ak-1)*sigma2e0;
		}
		}

	double ae   = 0.1;
	double be   = 0.1;
	double a_b  = 0.1;
	double b_b  = 0.1;
	//double a0   = 400;
	//double b0   = 40;
	double a0   = 1;
	double b0   = 0.1;
	double a_lambda = a0 + n_k;
	double b_lambda = b0;
	double Elogsigmae, tx_ywx_res, xtxabk, A;

	////initial values for a_k, b_k and sigma2k
	Ebeta2k = (mik_beta.cwiseProduct(mik_beta)) + sik2_beta;
	Ebeta   = (mik_beta.cwiseProduct(pik_beta)).rowwise().sum();	
	for (size_t k=0; k<n_k; k++)  {
		a_k(k) = (pik_beta.col(k).sum())/2 + ak;
		b_k(k) = (Ebeta2k.col(k).sum())*(a_e/b_e)/2 + bk(k);
		sigma2k(k) = b_k(k)/(a_k(k)-1);
		}

	A = (y_res).dot(y_res);
	B1=(1/sigma2k.array()).rowwise().replicate(n_snp);
	B1.col(0).fill(0);

	Ebeta2k = Ebeta2k.cwiseProduct(B1.transpose());
	double B = Ebeta2k.sum();
	a_e = n_idv + n_snp*0.1 + ae;
	b_e = (A + B + 2*be)/2;
	
	size_t int_step = 0;
	size_t max_step = n_idv*sqrt(10);
	double delta = 10;
	VectorXd ELBO = VectorXd::Zero(max_step);
	
	//WritelmmBeta(beta);
	if (n_k == 1)   {
		WritelmmBeta(beta);
		pheno_mean = m_alpha(0);
		}

	/////////////////// variational bayesian/////////////
	else {
	while ((int_step<max_step) && (delta > 1e-7)) {

		Elogsigmae = 0.5*(log(b_e)-gsl_sf_psi(a_e));
		//////////////  sampling the mixture snp effects 
		XEbeta.setZero(); //set Gbeta to zero first
			for (size_t i=0; i<n_snp; i++)   {
				x_col   =  UtX.col(i);
				XEbeta += x_col*Ebeta(i);
		}

		lambda_k(n_k-1) = 0;
		double ab_e = a_e/b_e; // E(sigma2e-2)
		double ba_e = b_e/a_e; // 1/E(sigma2e-2)

		y_res.setZero();
		y_res = Uty - WEalpha - Utu;
		for (size_t i=0; i<n_snp; i++)   {
			x_col      = UtX.col(i);
			XEbeta    -= x_col*Ebeta(i);
			tx_ywx_res = x_col.dot(y_res - XEbeta);
			for (size_t k=0; k<n_k; k++)   {
				if (k == 0) {
					mik_beta(i,k)  = 0; 
					sik2_beta(i,k) = 0; 
					pikexp1(k)     = 0;
					Elogsigmak(k)  = 0;
					}
				else {
					xtxabk         = xtx(i) + a_k(k)/b_k(k);
					mik_beta(i,k)  = tx_ywx_res/xtxabk;
					sik2_beta(i,k) = ba_e/xtxabk;
					Elogsigmak(k)  = 0.5*(log(b_k(k))-gsl_sf_psi(a_k(k)));
					pikexp1(k)     = pow(mik_beta(i,k),2)/(2*sik2_beta(i,k)) + log(sqrt(sik2_beta(i,k))) - Elogsigmae - Elogsigmak(k);
					//pikexp1(k)     =-pow(mik_beta(i,k),2)*(a_k(k)/b_k(k))*ab_e/2- Elogsigmae - Elogsigmak(k);
				}

				//if (k == (n_k-1))   {
					//vk(k)	  = 1;
					//Elogvk(k) = 0;
					//}
				//else   {
					Elogvk(k) = gsl_sf_psi(kappa_k(k))-gsl_sf_psi(kappa_k(k)+lambda_k(k));
				//}
			
				sumElogvl(k)  = sum_Elogvl(lambda_k,kappa_k,k);
				pikexp2(k)    = Elogvk(k) + sumElogvl(k);
				}

			index = pikexp1 + pikexp2;
			index = (index.array() - index.maxCoeff()).exp();
			pik_beta.row(i) = index/index.sum();


			/*
			/// sort sik2_beta and pik_beta
			gsl_vector *index0 = gsl_vector_alloc (n_k);			
			gsl_vector *sik2_beta0 = gsl_vector_alloc (n_k);
			for (size_t k=0; k<n_k; k++)   {
				gsl_vector_set(index0,k,index(k));
				gsl_vector_set(sik2_beta0,k,sik2_beta(i,k));
				}
			//gsl_sort_vector(index0);
			//gsl_vector_reverse(index0);
			gsl_permutation * perm0 = gsl_permutation_alloc (n_k);
			gsl_vector *perm1 = gsl_vector_alloc (n_k);
			VectorXd perm2(n_k);
			gsl_permutation_init (perm0);
			gsl_sort_vector_index (perm0, sik2_beta0);
			gsl_sort_vector(sik2_beta0);
			for (size_t k=0; k<n_k; k++)   {
				sik2_beta(i,k) = gsl_vector_get(sik2_beta0,k);
				gsl_vector_set(perm1, k, perm0->data[k]);
				perm2(k) = gsl_vector_get(perm1,k);
				index(k) = gsl_vector_get(index0,perm2(k));
				}*/

			//pik_beta.row(i) = index/index.sum();
			Ebeta(i) = ((mik_beta.row(i)).dot(pik_beta.row(i)));
			XEbeta  += x_col*Ebeta(i);
			//cout<<setprecision(7)<<sik2_beta.row(i)<<" ";
			//cout<<" "<<perm2(0)<<" "<<perm2(1)<<" "<<perm2(2)<<" "<<perm2(3)<<endl;
			//cout<<" "<<pik_beta.row(i)<<endl;
			/// sort sik2_beta and pik_beta
		}
		//////////////  sampling the mixture snp effects 

		WEalpha.setZero();
		for (size_t j=0;j<n_j; j++)   {
			WEalpha += UtW.col(j)*Ealpha(j);
		}
		
		y_res.setZero();
		y_res = Uty - XEbeta - Utu;
		for (size_t j=0;j<n_j; j++)   {
			WEalpha    -= UtW.col(j)*Ealpha(j);
			Ealpha(j)   = 1/wtw(j)*(UtW.col(j).dot(y_res - WEalpha));
			s2_alpha(j) = ba_e/wtw(j);
			WEalpha    += UtW.col(j)*Ealpha(j);
			}

		a_lambda = a0 + n_k;
		b_lambda = b0 - sum_b_lambda(lambda_k,kappa_k,n_k);

		for (size_t k=0; k<n_k-1; k++)   {
			 kappa_k(k) = pik_beta.col(k).sum() + 1;
			lambda_k(k) = sum_lambda_k(pik_beta,k,n_k) + a_lambda/b_lambda;
			}

		Ebeta2k = mik_beta.cwiseProduct(mik_beta) + sik2_beta;
		Ebeta2k =  Ebeta2k.cwiseProduct(pik_beta);

		for (size_t k=0; k<n_k; k++)   {
			a_k(k) = (pik_beta.col(k).sum())/2 + ak;
			b_k(k) =  (Ebeta2k.col(k).sum())*ab_e/2 + bk(k);
			}

		y_res.setZero();
		y_res = Uty - XEbeta - WEalpha;
		double ab = a_b/b_b;
		for (size_t i=0; i<n_idv; i++)   {
			if (D(i)  == 0)   {
				V(i)  = 0;
				Utu(i)= 0;
				Ue(i) = 0; 
				Ub(i) = 0;
				}
			else   {
				double abD = ab/D(i);
				V(i)  = ba_e/(abD + 1);
				Utu(i)= y_res(i)/(abD + 1);
				Ue(i) = (Utu(i)*Utu(i)+V(i))*abD;       //for sigma2e
				Ub(i) = (Utu(i)*Utu(i)+V(i))*ab_e/D(i); //for sigma2b
			}
			//bv(i) = y_res(i)/(D(i) + a_b/b_b); 
			}

		VectorXd VarBeta;
		VarBeta = Ebeta2k.rowwise().sum()-Ebeta.cwiseProduct(Ebeta);

		A =(y_res-Utu).dot(y_res-Utu)+wtw.dot(s2_alpha)+V.sum()+xtx.dot(VarBeta);
		B1=(a_k.array()/b_k.array()).rowwise().replicate(n_snp);
		B1.col(0).fill(0);
		double B  = (Ebeta2k.cwiseProduct(B1.transpose())).sum();
		double Gn = pik_beta.rightCols(n_k-1).sum();
		
		//mixture_no=Gn;

		a_e = n_idv + Gn/2 + ae;
		b_e = (A + B + Ue.sum())/2 + be;

		a_b = n_idv/2 + ae;
		b_b = Ub.sum()/2 + be;

		int_step++;
		VectorXd ab_k = a_k.array()/b_k.array();

		ELBO(int_step)=(lgamma(a_e)-a_e*log(b_e) +
						lgamma(a_b)-a_b*log(b_b) + a_b +
						ELBO1(a_k,b_k,n_k) +
						ELBO2(kappa_k,lambda_k,n_k) -
						ELBO3(pik_beta,sik2_beta) +
						0.5*( s2_alpha.array().log().sum()) + 
						0.5*((V.array()+1e-10).log().sum()) + 
						lgamma(a_lambda)-a_lambda*log(b_lambda))+a_lambda -
						be*(ab_k.tail(n_k-1).sum()+ab+a_lambda/b_lambda);

		////////////////////////////////////
		delta = abs((ELBO(int_step)-ELBO(int_step-1))/ELBO(int_step));

		//if ((int_step+1)%10==0) {cout<<int_step+1<<" "<<setprecision(5)<<delta<<" "<<ELBO(int_step)<<endl;}
		}
	}

	y_res.setZero();
	y_res = Uty - XEbeta - WEalpha;
	for (size_t i=0; i<n_idv; i++)   {
		bv(i) = y_res(i)/(D(i) + a_b/b_b); 
	}
	cout<<"Compute the polygenic effects ..."<<endl;
	cout<<"variational bayes is finished"<<endl;
	/*
	gsl_vector *gsl_bv    = gsl_vector_alloc (n_idv);
	gsl_vector *gsl_alpha = gsl_vector_alloc (n_snp);
	gsl_vector_set_zero (gsl_bv);
	VectorXd eigen_alpha = VectorXd::Zero (n_snp);
	for (size_t i=0; i<n_idv; i++) {
		gsl_vector_set (gsl_bv,i,bv(i));
	}
	gsl_blas_dgemv (CblasTrans,1.0/(double)n_snp,UtX,gsl_bv,0.0,gsl_alpha);
	for (size_t i = 0; i < n_snp; i++) {
		eigen_alpha(i) = gsl_vector_get(gsl_alpha, i);
		}
	*/

	VectorXd eigen_alpha = VectorXd::Zero (n_snp);
	eigen_alpha = UtX.transpose()*bv/n_snp;
	//WritevbBeta (Ebeta+eigen_alpha);
	WriteCoeff(eigen_alpha, Ebeta);
	pheno_mean = Ealpha(0);
	       LBO = ELBO(int_step);

	cout << "    delta = "<<scientific<<setprecision(6)<<delta << endl;
	cout << "iteraions = "<<scientific<<int_step+1<<endl;
	cout << "     ELBO = "<<scientific<<setprecision(6)<<ELBO(int_step) << endl;
	cout << "Computaion Time for DPR.VB = ";
	cout << setprecision(3)<<(clock()-time_begin)/(double(CLOCKS_PER_SEC)*60.0);
	cout <<" min"<<endl;
	cout << endl;

	////////////
	xtx.resize(0);xty.resize(0);wtw.resize(0);wty.resize(0);
	Utu.resize(0);Ue.resize(0);Ub.resize(0);bv.resize(0);
	V.resize(0);M.resize(0);Ealpha.resize(0);post_Ealpha.resize(0);
	m_alpha.resize(0);s2_alpha.resize(0);WEalpha.resize(0);
	XEbeta.resize(0);Ebeta.resize(0);post_Ebeta.resize(0);
	Ebeta2k.resize(0,0);B1.resize(0,0);mik_beta.resize(0,0);
	sik2_beta.resize(0,0);pik_beta.resize(0,0);beta_beta.resize(0,0);
	gamma_beta.resize(0,0);sigma2k.resize(0);vk.resize(0);
	Elogsigmak.resize(0);Elogvk.resize(0);sumElogvl.resize(0);
	pikexp1.resize(0);pikexp2.resize(0);index.resize(0);
	a_k.resize(0);b_k.resize(0);kappa_k.resize(0);
	lambda_k.resize(0);y_res.resize(0);bk.resize(0);
	eigen_alpha.resize(0);
}


void LDR::Gibbs (const MatrixXd UtX, const VectorXd Uty, const MatrixXd UtW,const VectorXd D, const VectorXd Wbeta, const VectorXd se_Wbeta, const VectorXd beta, double lambda) 
{
	
	setprecision(5); 
	clock_t time_begin  = clock();
	size_t n_snp = UtX.cols();
	size_t n_idv = Uty.size();
	size_t n_j   = UtW.cols();
	VectorXd x_col(n_idv);

	VectorXd xtx(n_snp); 
	VectorXd xty(n_snp); 
	VectorXd wtw(n_j);		 
	VectorXd wty(n_j);
	VectorXd Utu = VectorXd::Zero(n_idv);
	VectorXd Ue(n_idv);
	VectorXd Ub(n_idv);
	VectorXd bv0 = VectorXd::Zero(n_idv); //for computation of alpha
	VectorXd bv  = VectorXd::Zero(n_idv); //for computation of alpha
	VectorXd V = VectorXd::Zero (n_idv);//compute Utu 
	VectorXd M = VectorXd::Zero (n_idv);//compute Utu
	VectorXd snp_label(n_snp);

	VectorXd Ealpha      = Wbeta;//intercept
	VectorXd post_Ealpha = VectorXd::Zero(n_j);//save intercept
	VectorXd m_alpha     = Wbeta;
	VectorXd s2_alpha    = se_Wbeta.array()*se_Wbeta.array();

	VectorXd WEalpha(n_idv); 
	double   sigma2e=0, lambdax;
	VectorXd XEbeta(n_idv); //G*beta
	VectorXd Ebeta(n_snp);
	VectorXd post_Ebeta(n_snp);//save beta
	MatrixXd Ebeta2k;
	MatrixXd B1;
	MatrixXd mik_beta = MatrixXd::Zero(n_snp, n_k);

	for (size_t i=0;i<n_snp; i++)  {
		//Eigen_matrix_get_row (UtX, i, x_col);
		x_col  = UtX.col(i);
		xtx(i) = x_col.dot(x_col);
		xty(i) = x_col.dot(Uty);
	}

	for (size_t j=0;j<n_j; j++)  {
		wtw(j) = UtW.col(j).dot(UtW.col(j));
		wty(j) = UtW.col(j).dot(Uty);
	}

	//// initial values for SNP and each normal
	//// component has the same initial values
	for (size_t i=0;i<n_snp; i++)  {
		for (size_t k=1;k<n_k; k++)  {
			mik_beta(i,k) = beta(i);
		}
	}

	mik_beta.col(0).fill(0);
	MatrixXd sik2_beta  =(MatrixXd::Random(n_snp, n_k)).array().abs();
	MatrixXd pik_beta   = MatrixXd::Zero  (n_snp, n_k);
	pik_beta = pik_beta.array() + (double)1/n_k;

	MatrixXd beta_beta  = MatrixXd::Zero(n_snp, n_k);
	MatrixXd gamma_beta = MatrixXd::Zero(n_snp, n_k);
	MatrixXd post_gamma_beta = gamma_beta;
	gamma_beta = gamma_beta.array() + (double)1/n_k;

	VectorXd sigma2k (n_k);
	VectorXd vk (n_k);
	VectorXd Elogsigmak (n_k);
	VectorXd Elogvk (n_k);
	VectorXd sumElogvl (n_k);

	VectorXd pikexp1(n_k);
	VectorXd pikexp2(n_k);
	VectorXd index;
	VectorXd a_k = (VectorXd::Random(n_k)).array().abs();
	VectorXd b_k = (VectorXd::Random(n_k)).array().abs();
	VectorXd  kappa_k = pik_beta.colwise().sum();
	VectorXd lambda_k = (VectorXd::Random(n_k)).array().abs();

	XEbeta.setZero();//set to zero
	for (size_t i=0; i<n_snp; i++)   {
		x_col   =  UtX.col(i);
		XEbeta += x_col*beta(i);
		}

	WEalpha.setZero(); //set to zero
	for (size_t j=0;j<n_j; j++)   {
		WEalpha += UtW.col(j)*m_alpha(j);
		}

	VectorXd y_res  = Uty - WEalpha - XEbeta;
	double sigma2e0 = ((y_res).dot(y_res))/(n_idv-1);
	cout << endl << endl;
	double a_e  = (2*n_idv + n_snp)/2;
	double b_e  = (a_e-1)*sigma2e0;

	double ak   = 21;
	if (lambda < 0.01)     {lambda =    0.01;}
	else if (lambda > 100) {lambda =     100;}
	else                   {lambda =  lambda;}
	double bk0  = lambda*(ak-1)/n_snp;

	VectorXd bk(n_k);
	bk(0) = bk0;
	for (size_t i=1;i<n_k; i++)  {
		bk(i) = bk(i-1)*1.7*sqrt(pow(i,i));
		}

	for (size_t i=0;i<n_snp; i++)  {
		for (size_t k=0;k<n_k; k++)  {
			sik2_beta(i,k) = bk(k)/(ak-1)*sigma2e0;
		}
	}

	double ae   = 0.1;
	double be   = 0.1;
	double a0   = 400;
	double b0   = 40;
	double a_lambda = a0 + n_k;
	double b_lambda = b0;
	double Elogsigmae, tx_ywx_res, xtxabk, A, post_Gn=0;
	double sigma2b = 0.001;

	////initial values for a_k, b_k and sigma2k
	Ebeta2k = (mik_beta.cwiseProduct(mik_beta)) + sik2_beta;
	for (size_t k=0; k<n_k; k++)  {
		a_k(k) = (pik_beta.col(k).sum())/2 + ak;
		b_k(k) = (Ebeta2k.col(k).sum())*(a_e/b_e)/2 + bk(k);
		sigma2k(k) = b_k(k)/(a_k(k)-1);
		}

	A = (y_res).dot(y_res);
	B1=(1/sigma2k.array()).rowwise().replicate(n_snp);
	B1.col(0).fill(0);

	Ebeta2k = Ebeta2k.cwiseProduct(B1.transpose());
	double B = Ebeta2k.sum();
	a_e = n_idv + n_snp*0.1 + ae;
	b_e = (A + B + 2*be)/2;
	////random seed
	gsl_rng * rs = gsl_rng_alloc (gsl_rng_mt19937);

	if (n_k == 1)   {
		WritelmmBeta(beta);
		pheno_mean = m_alpha(0);
		}

	// //  begin MCMC sampling
	else {
	for (size_t S=0; S<(w_step + s_step); S++)   { 

		sigma2e    = 1/gsl_ran_gamma(rs,a_e,1/b_e);
		Elogsigmae = log(sqrt(sigma2e));

		// save Ebeta for the mixture normal component
		Ebeta = (beta_beta.cwiseProduct(gamma_beta)).rowwise().sum();		
		
		if (S > (w_step-1))   {post_Ebeta += Ebeta;}

		//sample sigma2k and compute related quantities
		for (size_t k=0; k<n_k; k++)   {
			if (k==0) {
				sigma2k(k)    = 0; 
				Elogsigmak(k) = 0;
				}

			else {
				sigma2k(k)    = 1/gsl_ran_gamma(rs,a_k(k),1/b_k(k));
				Elogsigmak(k) = log(sqrt(sigma2k(k)));
				}

				if (k == (n_k-1)) {vk(k) = 0; Elogvk(k) = 0;}
				else {
					vk(k)     = gsl_ran_beta(rs,kappa_k(k),lambda_k(k));
					Elogvk(k) = log(vk(k));
				}
				sumElogvl(k)  = sum_Elogvl2(vk,k);
			}

		//////////////  sampling the mixture snp effects 
		XEbeta.setZero(); //set Gbeta to zero first
			for (size_t i=0; i<n_snp; i++)   {
				x_col   =  UtX.col(i);
				XEbeta += x_col*Ebeta(i);
		}

		y_res.setZero();
		y_res = Uty - WEalpha - Utu;

		if (S<(w_step + w_step)*tp)   {  

		//////////////////////////////////////
		for (size_t i=0; i<n_snp; i++)   {
			x_col      = UtX.col(i);
			XEbeta    -= x_col*Ebeta(i);
			tx_ywx_res = x_col.dot(y_res - XEbeta);
			for (size_t k=0; k<n_k; k++)   {

				if (k == 0) {
					mik_beta(i,k)  = 0; 
					sik2_beta(i,k) = 0; 
					pikexp1(k)     = 0;
					}

				else {
					xtxabk         = xtx(i) + 1/sigma2k(k);
					mik_beta(i,k)  = tx_ywx_res/xtxabk;
					sik2_beta(i,k) = sigma2e/xtxabk;
					pikexp1(k)     = mik_beta(i,k)*mik_beta(i,k)/(2*sik2_beta(i,k)) + log(sqrt(sik2_beta(i,k))) - Elogsigmae - Elogsigmak(k);
				}

				beta_beta(i,k) = gsl_ran_gaussian(rs,sqrt(sik2_beta(i,k))) + mik_beta(i,k);
				pikexp2(k)     = Elogvk(k) + sumElogvl(k);
				}

			index = pikexp1 + pikexp2;
			index = (index.array() - index.maxCoeff()).exp();
			pik_beta.row(i) = index/index.sum();

			/// sort sik2_beta and pik_beta
			/*gsl_vector *index0 = gsl_vector_alloc (n_k);			
			gsl_vector *sik2_beta0 = gsl_vector_alloc (n_k);
			for (size_t k=0; k<n_k; k++)   {
				gsl_vector_set(index0,k,index(k));
				gsl_vector_set(sik2_beta0,k,sik2_beta(i,k));
				}
			//gsl_sort_vector(index0);
			//gsl_vector_reverse(index0);
			gsl_permutation * perm0 = gsl_permutation_alloc (n_k);
			gsl_vector *perm1 = gsl_vector_alloc (n_k);
			VectorXd perm2(n_k);
			gsl_permutation_init (perm0);
			gsl_sort_vector_index (perm0, sik2_beta0);
			gsl_sort_vector(sik2_beta0);
			for (size_t k=0; k<n_k; k++)   {
				sik2_beta(i,k) = gsl_vector_get(sik2_beta0,k);
				gsl_vector_set(perm1, k, perm0->data[k]);
				perm2(k) = gsl_vector_get(perm1,k);
				index(k) = gsl_vector_get(index0,perm2(k));
				}
			pik_beta.row(i) = index/index.sum();
			//cout<<setprecision(7)<<sik2_beta.row(i)<<" ";
			//cout<<" "<<perm2(0)<<" "<<perm2(1)<<" "<<perm2(2)<<" "<<perm2(3)<<endl;
			//cout<<" "<<pik_beta.row(i)<<endl;
			/// sort sik2_beta and pik_beta			/// sort sik2_beta and pik_beta
*/
		// multinomial sampling
		double mult_prob[n_k];
		unsigned int mult_no[n_k];
		for (size_t k=0; k<n_k; k++)   {
			mult_prob[k] = pik_beta(i,k);
			}
		gsl_ran_multinomial(rs, n_k, 1, mult_prob, mult_no);
		for (size_t k=0; k<n_k; k++)   {
			gamma_beta(i,k) = mult_no[k];
		}
		Ebeta(i) = ((beta_beta.row(i)).dot(gamma_beta.row(i)));
		XEbeta  += x_col*Ebeta(i);
		}
		}
		//////////////////////////////////////
		
		else  {

		snp_label = permute(n_snp);
		size_t i, non_zero = 0, j = 0, T0 = n_snp*sp+1;
		////////////////////////////////////// 
		while (j<T0)  {
			i = snp_label(j);
			x_col      =  UtX.col(i);
			XEbeta    -= x_col*Ebeta(i);
			tx_ywx_res = x_col.dot(y_res - XEbeta);
			for (size_t k=0; k<n_k; k++)   {				
				if (k == 0)   {
					mik_beta(i,k)  = 0; 
					sik2_beta(i,k) = 0; 
					pikexp1(k)     = 0;
					}
				else {
					xtxabk         = xtx(i) + 1/sigma2k(k);
					mik_beta(i,k)  = tx_ywx_res/xtxabk;
					sik2_beta(i,k) = sigma2e/xtxabk;
					pikexp1(k)     = mik_beta(i,k)*mik_beta(i,k)/(2*sik2_beta(i,k)) + log(sqrt(sik2_beta(i,k))) - Elogsigmae - Elogsigmak(k);
				}
				beta_beta(i,k) = gsl_ran_gaussian(rs,sqrt(sik2_beta(i,k))) + mik_beta(i,k);
				pikexp2(k)     = Elogvk(k) + sumElogvl(k);
				}
			index = pikexp1 + pikexp2;
			index = (index.array() - index.maxCoeff()).exp();
			pik_beta.row(i) = index/index.sum();

		// multinomial sampling
		double mult_prob[n_k];
		unsigned int mult_no[n_k];
		for (size_t k=0; k<n_k; k++)   {
			mult_prob[k] = pik_beta(i,k);
			}
		gsl_ran_multinomial(rs, n_k, 1, mult_prob, mult_no);
		for (size_t k=0; k<n_k; k++)   {
			gamma_beta(i,k) = mult_no[k];
		}
		// multinomial sampling

		non_zero+= gamma_beta.row(i).tail(n_k-1).sum();
		Ebeta(i) = ((beta_beta.row(i)).dot(gamma_beta.row(i)));
		XEbeta  += x_col*Ebeta(i);
		j++;
		}
		}
		//////////////////////////////////////

		if (S > (w_step-1))   {post_gamma_beta += gamma_beta;}
		//////////////  sampling the mixture snp effects 

		WEalpha.setZero();
		for (size_t j=0;j<n_j; j++)   {
			WEalpha += UtW.col(j)*Ealpha(j);
		}
		
		y_res.setZero();
		y_res = Uty - XEbeta - Utu;
		for (size_t j=0;j<n_j; j++)   {
			WEalpha    -= UtW.col(j)*Ealpha(j);
			m_alpha(j)  = 1/wtw(j)*(UtW.col(j).dot(y_res - WEalpha));
			s2_alpha(j) = sigma2e/wtw(j);
			Ealpha(j)   = m_alpha(j) + gsl_ran_gaussian(rs,sqrt(s2_alpha(j)));
			WEalpha    += UtW.col(j)*Ealpha(j);
			}

		if (S>(w_step-1))   {post_Ealpha += Ealpha;}

		a_lambda = a0 + n_k;
		b_lambda = b0 - sum_Elogvl2(vk,n_k-1);
		lambdax  = gsl_ran_gamma(rs,a_lambda,1/b_lambda);

		for (size_t k=0; k<n_k-1; k++)   {
			kappa_k(k)  = gamma_beta.col(k).sum() + 1;
			lambda_k(k) = sum_lambda_k(gamma_beta,k,n_k) + lambdax;
			}

		Ebeta2k = beta_beta.cwiseProduct(beta_beta).cwiseProduct(gamma_beta);
		for (size_t k=0; k<n_k; k++)   {
			a_k(k) = (gamma_beta.col(k).sum())/2 + ak;
			b_k(k) = (Ebeta2k.col(k).sum())/(2*sigma2e) + bk(k);
			}
		y_res.setZero();
		y_res = Uty - XEbeta - WEalpha;
		for (size_t i=0; i<n_idv; i++) {
			V(i)  = sigma2b*D(i)/(sigma2b*D(i) + 1);
			M(i)  = y_res(i)*V(i);
			Utu(i)= M(i) + gsl_ran_gaussian(rs,sqrt(V(i)*sigma2e));
			if (D(i)  == 0) {
				Ue(i) = 0; 
				Ub(i) = 0;
				}
			else {
				Ue(i) = Utu(i)/(sigma2b*D(i)); //for sigma2e
				Ub(i) = Utu(i)/(sigma2e*D(i)); //for sigma2b
			}
			bv0(i) = y_res(i)*sigma2b/(sigma2b*D(i) + 1); 
			}

		if (S>(w_step-1))   {bv += bv0;}

		A  = (y_res-Utu).dot(y_res-Utu);
		B1 = (1/sigma2k.tail(n_k-1).array()).rowwise().replicate(n_snp);
		Ebeta2k   = Ebeta2k.rightCols(n_k-1).cwiseProduct(B1.transpose());
		double B  = Ebeta2k.sum();
		double Gn = gamma_beta.rightCols(n_k-1).sum();
		
		a_e = n_idv + Gn/2 + ae;
		b_e = (A + B + Ue.dot(Utu) + 2*be)/2;
		sigma2b = 1/gsl_ran_gamma(rs,n_idv/2+ae,1/(2*Ub.dot(Utu))+be);

		if (S > (w_step-1))   {post_Gn += Gn;}
		////////////////////////////////////
		//cout<<setprecision(6)<<sigma2b<<" "<<sigma2e<<" "<<sigma2b*sigma2e<<endl;

		size_t xstep = (s_step + w_step)*0.01+1;
		if (S%xstep==0||S==(s_step + w_step-1)) {ProgressBar ("MCMC sampling ", S, s_step+w_step-1);}
		}
	}

	mixture_no = post_Gn/s_step;
	cout<<endl<<endl<<"MCMC sampling is finished and now compute the polygenic effects ..."<<endl;

	VectorXd eigen_alpha = VectorXd::Zero (n_snp);
	eigen_alpha = UtX.transpose()*(bv/s_step)/n_snp;
	WriteCoeff (eigen_alpha, post_Ebeta/s_step);
	pheno_mean = post_Ealpha(0)/s_step;

	cout<<endl<< "Computaion Time for LDR.Gibbs = ";
	cout<<(clock()-time_begin)/(double(CLOCKS_PER_SEC)*60.0)<<" min"<<endl<<endl;
/*
	//output mik and pik
	ofstream outfile1;
	ofstream outfile2;
	ofstream outfile3;
	ofstream outfile4;
	outfile1.open ("mik.txt", ios::out | ios::binary);
	outfile2.open ("pik.txt", ios::out | ios::binary);
	outfile3.open ("sik.txt", ios::out | ios::binary);
	//outfile4.open ("sigmak.txt", ios::out | ios::binary);
	
	for (size_t i = 0; i<n_snp; i++)  {
		for (size_t k = 0; k<n_k; k++)  {
			outfile1 << setprecision(4) <<  mik_beta(i,k)<<"\t";
			outfile2 << setprecision(4) <<  pik_beta(i,k)<<"\t";
			outfile3 << setprecision(4) << sik2_beta(i,k)<<"\t";
			//outfile4 << setprecision(4) << b_k(k)/a_k(k) <<"\t";
			}
		outfile1 << endl;
		outfile2 << endl;
		outfile3 << endl;
		//outfile4 << endl;
	}
	outfile1.close();
	outfile2.close();
	outfile3.close();
	//outfile4.close();
	//output mik and pik
*/

	////////////
	xtx.resize(0);xty.resize(0);wtw.resize(0);wty.resize(0);
	Utu.resize(0);Ue.resize(0);Ub.resize(0);bv0.resize(0);
	bv.resize(0);V.resize(0);M.resize(0);Ealpha.resize(0);
	post_Ealpha.resize(0);m_alpha.resize(0);s2_alpha.resize(0);
	WEalpha.resize(0);XEbeta.resize(0);Ebeta.resize(0);
	post_Ebeta.resize(0);Ebeta2k.resize(0,0);
	B1.resize(0,0);mik_beta.resize(0,0);sik2_beta.resize(0,0);
	pik_beta.resize(0,0);beta_beta.resize(0,0);gamma_beta.resize(0,0);
	sigma2k.resize(0);vk.resize(0);Elogsigmak.resize(0);
	Elogvk.resize(0);sumElogvl.resize(0);pikexp1.resize(0);
	pikexp2.resize(0);index.resize(0);a_k.resize(0);
	b_k.resize(0);kappa_k.resize(0);lambda_k.resize(0);
	y_res.resize(0);bk.resize(0);
	eigen_alpha.resize(0);
}


void LDR::Gibbs_without_u (const MatrixXd UtX, const VectorXd Uty, const MatrixXd UtW,const VectorXd D, const VectorXd Wbeta, const VectorXd se_Wbeta, const VectorXd beta, double lambda) 
{
	
	setprecision(5); 
	//clock_t time_begin  = clock();
	size_t n_snp = UtX.cols();
	size_t n_idv = Uty.size();
	size_t n_j   = UtW.cols();
	VectorXd x_col(n_idv);

	VectorXd xtx(n_snp);
	VectorXd xty(n_snp);
	VectorXd wtw(n_j);
	VectorXd wty(n_j);
	VectorXd Utu = VectorXd::Zero(n_idv);
	VectorXd Ue(n_idv);
	VectorXd Ub(n_idv);

	VectorXd bv0 = VectorXd::Zero(n_idv); //for computation of alpha
	VectorXd bv  = VectorXd::Zero(n_idv); //for computation of alpha
	VectorXd V = VectorXd::Zero (n_idv);//compute Utu 
	VectorXd M = VectorXd::Zero (n_idv);//compute Utu
	VectorXd snp_label(n_snp);

	VectorXd Ealpha      = Wbeta;//intercept
	VectorXd post_Ealpha = VectorXd::Zero(n_j);//save intercept
	VectorXd m_alpha     = Wbeta;
	VectorXd s2_alpha    = se_Wbeta.array()*se_Wbeta.array();

	VectorXd WEalpha(n_idv); 
	double   sigma2e=0, lambdax;
	VectorXd XEbeta(n_idv); //G*beta
	VectorXd Ebeta(n_snp);
	VectorXd post_Ebeta = VectorXd::Zero(n_snp);//save beta
	MatrixXd Ebeta2k;
	MatrixXd B1;
	MatrixXd mik_beta = MatrixXd::Zero(n_snp, n_k);

	VectorXd    llike = VectorXd::Zero(w_step+s_step);

	for (size_t i=0;i<n_snp; i++)  {
		//Eigen_matrix_get_row (UtX, i, x_col);
		x_col  = UtX.col(i);
		xtx(i) = x_col.dot(x_col);
		xty(i) = x_col.dot(Uty);
	}

	for (size_t j=0;j<n_j; j++)  {
		wtw(j) = UtW.col(j).dot(UtW.col(j));
		wty(j) = UtW.col(j).dot(Uty);
	}

	//// initial values for SNP and each normal
	//// component has the same initial values
	for (size_t i=0;i<n_snp; i++)  {
		for (size_t k=1;k<n_k; k++)  {
			mik_beta(i,k) = beta(i);
		}
	}

	mik_beta.col(0).fill(0);
	MatrixXd sik2_beta  =(MatrixXd::Random(n_snp, n_k)).array().abs();
	MatrixXd pik_beta   = MatrixXd::Zero  (n_snp, n_k);
	pik_beta = pik_beta.array() + (double)1/n_k;

	MatrixXd beta_beta  = MatrixXd::Zero(n_snp, n_k);
	MatrixXd gamma_beta = MatrixXd::Zero(n_snp, n_k);
	MatrixXd post_gamma_beta = gamma_beta;
	gamma_beta = gamma_beta.array() + (double)1/n_k;

	VectorXd sigma2k (n_k);
	VectorXd vk (n_k);
	VectorXd Elogsigmak (n_k);
	VectorXd Elogvk (n_k);
	VectorXd sumElogvl (n_k);

	VectorXd pikexp1(n_k);
	VectorXd pikexp2(n_k);
	VectorXd index;
	VectorXd a_k = (VectorXd::Random(n_k)).array().abs();
	VectorXd b_k = (VectorXd::Random(n_k)).array().abs();
	VectorXd  kappa_k = pik_beta.colwise().sum();
	VectorXd lambda_k = (VectorXd::Random(n_k)).array().abs();

	XEbeta.setZero();//set to zero
	for (size_t i=0; i<n_snp; i++)   {
		x_col   =  UtX.col(i);
		XEbeta += x_col*beta(i);
		}

	WEalpha.setZero(); //set to zero
	for (size_t j=0;j<n_j; j++)   {
		WEalpha += UtW.col(j)*m_alpha(j);
		}

	VectorXd y_res  = Uty - WEalpha - XEbeta;
	double sigma2e0 = ((y_res).dot(y_res))/(n_idv-1);
	cout << endl << endl;
	double a_e  = (2*n_idv + n_snp)/2;
	double b_e  = (a_e-1)*sigma2e0;

	double ak   = 21;
	if (lambda < 0.01)     {lambda =    0.01;}
	else if (lambda > 100) {lambda =     100;}
	else                   {lambda =  lambda;}
	double bk0  = lambda*(ak-1)/n_snp;

	VectorXd bk(n_k);
	bk(0) = bk0;
	for (size_t i=1;i<n_k; i++)  {
		bk(i) = bk(i-1)*1.7*sqrt(pow(i,i));
		}

	for (size_t i=0;i<n_snp; i++)  {
		for (size_t k=0;k<n_k; k++)  {
			sik2_beta(i,k) = bk(k)/(ak-1)*sigma2e0;
		}
	}

	double ae   = 0.1;
	double be   = 0.1;
	double a0   = 400;
	double b0   = 40;
	double a_lambda = a0 + n_k;
	double b_lambda = b0;
	double Elogsigmae, tx_ywx_res, xtxabk, A, post_Gn=0;
	double sigma2b = 0.1;

	double post_llike=0, post_sigma2e=0, post_sigma2b=0;

	VectorXd sigma2bX = (VectorXd::Random(w_step+s_step+1)).array().abs();
	sigma2bX(0)=sigma2b;

	VectorXd h = (VectorXd::Random(w_step+s_step+1)).array().abs();
	VectorXd H0(n_idv);
	VectorXd H(n_idv);

	//for (size_t k=0; k<(w_step+s_step+1); k++)  {
		//cout<<h(k)<<endl;
	//}

	////initial values for a_k, b_k and sigma2k
	Ebeta2k = (mik_beta.cwiseProduct(mik_beta)) + sik2_beta;
	for (size_t k=0; k<n_k; k++)  {
		a_k(k) = (pik_beta.col(k).sum())/2 + ak;
		b_k(k) = (Ebeta2k.col(k).sum())*(a_e/b_e)/2 + bk(k);
		sigma2k(k) = b_k(k)/(a_k(k)-1);
		}

	A = (y_res).dot(y_res);
	B1=(1/sigma2k.array()).rowwise().replicate(n_snp);
	B1.col(0).fill(0);

	Ebeta2k = Ebeta2k.cwiseProduct(B1.transpose());
	double B = Ebeta2k.sum();
	a_e = n_idv + n_snp*0.1 + ae;
	b_e = (A + B + 2*be)/2;

	////random seed
	gsl_rng * rs = gsl_rng_alloc (gsl_rng_mt19937);

	if (n_k == 1)   {
		WritelmmBeta(beta);
		pheno_mean = m_alpha(0);
		}

	// //  begin MCMC sampling
	else {
	for (size_t S=0; S<(w_step + s_step); S++)   {

		sigma2e    = 1/gsl_ran_gamma(rs,a_e,1/b_e);
		Elogsigmae = log(sqrt(sigma2e));

		// save Ebeta for the mixture normal component
		Ebeta = (beta_beta.cwiseProduct(gamma_beta)).rowwise().sum();

		if (S > (w_step-1))   {post_Ebeta += Ebeta;}

		//sample sigma2k and compute related quantities
		for (size_t k=0; k<n_k; k++)   {
			if (k==0) {
				sigma2k(k)    = 0; 
				Elogsigmak(k) = 0;
				}

			else {
				sigma2k(k)    = 1/gsl_ran_gamma(rs,a_k(k),1/b_k(k));
				Elogsigmak(k) = log(sqrt(sigma2k(k)));
				}

				if (k == (n_k-1)) {vk(k) = 0; Elogvk(k) = 0;}
				else {
					vk(k)     = gsl_ran_beta(rs,kappa_k(k),lambda_k(k));
					Elogvk(k) = log(vk(k));
				}
				sumElogvl(k)  = sum_Elogvl2(vk,k);
			}

		//////////////  sampling the mixture snp effects 
		XEbeta.setZero(); //set Gbeta to zero first
			for (size_t i=0; i<n_snp; i++)   {
				x_col   =  UtX.col(i);
				XEbeta += x_col*Ebeta(i);
		}

		y_res.setZero();
		y_res = Uty - WEalpha;

		H0 = (sigma2b*D).array()+1;
		H0 = 1/H0.array();

		if (S<(w_step + w_step)*tp)   {  

		//////////////////////////////////////
		for (size_t i=0; i<n_snp; i++)   {
			x_col      = UtX.col(i);
			XEbeta    -= x_col*Ebeta(i);

			H = x_col.array()*H0.array();
			tx_ywx_res = H.dot(y_res - XEbeta);

			for (size_t k=0; k<n_k; k++)   {

				if (k == 0) {
					mik_beta(i,k)  = 0; 
					sik2_beta(i,k) = 0; 
					pikexp1(k)     = 0;
					}

				else {
					xtxabk         = H.dot(x_col) + 1/sigma2k(k);
					mik_beta(i,k)  = tx_ywx_res/xtxabk;
					sik2_beta(i,k) = sigma2e/xtxabk;
					pikexp1(k)     = mik_beta(i,k)*mik_beta(i,k)/(2*sik2_beta(i,k)) + log(sqrt(sik2_beta(i,k))) - Elogsigmae - Elogsigmak(k);
				}

				beta_beta(i,k) = gsl_ran_gaussian(rs,sqrt(sik2_beta(i,k))) + mik_beta(i,k);
				pikexp2(k)     = Elogvk(k) + sumElogvl(k);
				}

			index = pikexp1 + pikexp2;
			index = (index.array() - index.maxCoeff()).exp();
			pik_beta.row(i) = index/index.sum();

			/// sort sik2_beta and pik_beta
			/*
			gsl_vector *index0 = gsl_vector_alloc (n_k);			
			gsl_vector *sik2_beta0 = gsl_vector_alloc (n_k);
			for (size_t k=0; k<n_k; k++)   {
				gsl_vector_set(index0,k,index(k));
				gsl_vector_set(sik2_beta0,k,sik2_beta(i,k));
				}
			//gsl_sort_vector(index0);
			//gsl_vector_reverse(index0);
			gsl_permutation * perm0 = gsl_permutation_alloc (n_k);
			gsl_vector *perm1 = gsl_vector_alloc (n_k);
			VectorXd perm2(n_k);
			gsl_permutation_init (perm0);
			gsl_sort_vector_index (perm0, sik2_beta0);
			gsl_sort_vector(sik2_beta0);
			for (size_t k=0; k<n_k; k++)   {
				sik2_beta(i,k) = gsl_vector_get(sik2_beta0,k);
				gsl_vector_set(perm1, k, perm0->data[k]);
				perm2(k) = gsl_vector_get(perm1,k);
				index(k) = gsl_vector_get(index0,perm2(k));
				}
			pik_beta.row(i) = index/index.sum();
			//cout<<setprecision(7)<<sik2_beta.row(i)<<" ";
			//cout<<" "<<perm2(0)<<" "<<perm2(1)<<" "<<perm2(2)<<" "<<perm2(3)<<endl;
			//cout<<" "<<pik_beta.row(i)<<endl;
			/// sort sik2_beta and pik_beta			/// sort sik2_beta and pik_beta
*/

		// multinomial sampling
		double mult_prob[n_k];
		unsigned int mult_no[n_k];

		for (size_t k=0; k<n_k; k++)   {
			mult_prob[k] = pik_beta(i,k);
			}

		gsl_ran_multinomial(rs, n_k, 1, mult_prob, mult_no);
		for (size_t k=0; k<n_k; k++)   {
			gamma_beta(i,k) = mult_no[k];
		}

		Ebeta(i) = ((beta_beta.row(i)).dot(gamma_beta.row(i)));
		XEbeta  += x_col*Ebeta(i);
		}

		}
		//////////////////////////////////////
		
		else  {

		snp_label = permute(n_snp);
		size_t i, non_zero = 0, j = 0, T0 = n_snp*sp+1;
		////////////////////////////////////// 
		while (j<T0)  {
			i = snp_label(j);
			x_col      =  UtX.col(i);
			XEbeta    -= x_col*Ebeta(i);
			tx_ywx_res = x_col.dot(y_res - XEbeta);
			for (size_t k=0; k<n_k; k++)   {				
				if (k == 0)   {
					mik_beta(i,k)  = 0; 
					sik2_beta(i,k) = 0; 
					pikexp1(k)     = 0;
					}
				else {
					xtxabk         = xtx(i) + 1/sigma2k(k);
					mik_beta(i,k)  = tx_ywx_res/xtxabk;
					sik2_beta(i,k) = sigma2e/xtxabk;
					pikexp1(k)     = mik_beta(i,k)*mik_beta(i,k)/(2*sik2_beta(i,k)) + log(sqrt(sik2_beta(i,k))) - Elogsigmae - Elogsigmak(k);
				}
				beta_beta(i,k) = gsl_ran_gaussian(rs,sqrt(sik2_beta(i,k))) + mik_beta(i,k);
				pikexp2(k)     = Elogvk(k) + sumElogvl(k);
				}
			index = pikexp1 + pikexp2;
			index = (index.array() - index.maxCoeff()).exp();
			pik_beta.row(i) = index/index.sum();

		// multinomial sampling
		double mult_prob[n_k];
		unsigned int mult_no[n_k];
		for (size_t k=0; k<n_k; k++)   {
			mult_prob[k] = pik_beta(i,k);
			}
		gsl_ran_multinomial(rs, n_k, 1, mult_prob, mult_no);
		for (size_t k=0; k<n_k; k++)   {
			gamma_beta(i,k) = mult_no[k];
		}
		// multinomial sampling

		non_zero+= gamma_beta.row(i).tail(n_k-1).sum();
		Ebeta(i) = ((beta_beta.row(i)).dot(gamma_beta.row(i)));
		XEbeta  += x_col*Ebeta(i);
		j++;
		}
		}
		//////////////////////////////////////

		if (S > (w_step-1))   {post_gamma_beta += gamma_beta;}
		//////////////  sampling the mixture snp effects 

		WEalpha.setZero();
		for (size_t j=0;j<n_j; j++)   {
			WEalpha += UtW.col(j)*Ealpha(j);
		}

		y_res.setZero();
		y_res = Uty - XEbeta;
		for (size_t j=0;j<n_j; j++)   {
			WEalpha    -= UtW.col(j)*Ealpha(j);

					  H = UtW.col(j).array()*H0.array();

			m_alpha(j)  = (H.dot(y_res - WEalpha))/(H.dot(UtW.col(j)));
			s2_alpha(j) = sigma2e/(H.dot(UtW.col(j)));
			Ealpha(j)   = m_alpha(j) + gsl_ran_gaussian(rs,sqrt(s2_alpha(j)));
			WEalpha    += UtW.col(j)*Ealpha(j);
			}

		if (S>(w_step-1))   {post_Ealpha += Ealpha;}

		a_lambda = a0 + n_k;
		b_lambda = b0 - sum_Elogvl2(vk,n_k-1);
		lambdax  = gsl_ran_gamma(rs,a_lambda,1/b_lambda);

		for (size_t k=0; k<n_k-1; k++)   {
			kappa_k(k)  = gamma_beta.col(k).sum() + 1;
			lambda_k(k) = sum_lambda_k(gamma_beta,k,n_k) + lambdax;
			}

		Ebeta2k = beta_beta.cwiseProduct(beta_beta).cwiseProduct(gamma_beta);
		for (size_t k=0; k<n_k; k++)   {
			a_k(k) = (gamma_beta.col(k).sum())/2 + ak;
			b_k(k) = (Ebeta2k.col(k).sum())/(2*sigma2e) + bk(k);
			}

		/*y_res.setZero();
		y_res = Uty - XEbeta - WEalpha;
		for (size_t i=0; i<n_idv; i++) {
			V(i)  = sigma2b*D(i)/(sigma2b*D(i) + 1);
			M(i)  = y_res(i)*V(i);
			Utu(i)= M(i) + gsl_ran_gaussian(rs,sqrt(V(i)*sigma2e));
			if (D(i)  == 0) {
				Ue(i) = 0; 
				Ub(i) = 0;
				}
			else {
				Ue(i) = Utu(i)/(sigma2b*D(i)); //for sigma2e
				Ub(i) = Utu(i)/(sigma2e*D(i)); //for sigma2b
			}
			bv0(i) = y_res(i)*sigma2b/(sigma2b*D(i) + 1); 
			}

		if (S>(w_step-1))   {bv += bv0;}
		*/

		y_res.setZero();
		y_res = Uty - XEbeta - WEalpha;

		bv0 = (y_res*sigma2b).array()*H0.array();
		if (S>(w_step-1))   {bv += bv0;}

		H = y_res.array()*H0.array();

		A  = H.dot(y_res);
		B1 = (1/sigma2k.tail(n_k-1).array()).rowwise().replicate(n_snp);
		Ebeta2k   = Ebeta2k.rightCols(n_k-1).cwiseProduct(B1.transpose());
		double B  = Ebeta2k.sum();
		double Gn = gamma_beta.rightCols(n_k-1).sum();
		
		a_e = n_idv/2 + Gn/2 + ae;
		b_e = (A + B + 2*be)/2;

		/*double sigma2bX_new = sigma2bX(S) + gsl_ran_gaussian(rs,sqrt(1));
		double ratio1 = log_sigma2b(D,y_res,sigma2bX_new,sigma2e,ae,be);
		double ratio2 = log_sigma2b(D,y_res,sigma2bX(S), sigma2e,ae,be);
		double ratio  = exp(ratio1-ratio2);
		if (ratio>1) {ratio = 1;}
		else         {ratio = ratio;}
		if (sigma2bX_new<0) {ratio = 0;}
		*/

		//double h_new  = gsl_rng_uniform(rs);
		double h_new  = gsl_ran_beta(rs,2,8);
		double ratio1 = log_h(D,y_res,h_new,sigma2e,ae,be)-log(gsl_ran_beta_pdf(h_new,2,8));
		double ratio2 = log_h(D,y_res,h(S), sigma2e,ae,be)-log(gsl_ran_beta_pdf(h(S),2,8));
		double ratio  = exp(ratio1-ratio2);
		if (ratio>1) {ratio = 1;}
		else         {ratio = ratio;}
		double u = gsl_rng_uniform(rs);
		if (u<ratio) {h(S+1) = h_new;}
		else         {h(S+1) = h(S);}

		sigma2b = h(S+1)/(1-h(S+1));

		//cout<<S<<" "<<setprecision(6)<<h_new<<" "<<sigma2b<<" "<<sigma2e<<" "<<sigma2b*sigma2e<<endl;
		//cout<<setprecision(6)<<sigma2bX_new<<" "<<sigma2b<<" "<<sigma2e<<" "<<sigma2b*sigma2e<<endl;
		if (S > (w_step-1))   {post_Gn += Gn;}
		////////////////////////////////////

		double llike0=logLike(D,y_res,sigma2b,sigma2e); 
		     llike(S)=llike0; 

		if (S > (w_step-1))   {	
			post_llike   += llike0;
			post_sigma2e += sigma2e;
			post_sigma2b += sigma2b;
			}

		size_t xstep = (s_step + w_step)*0.2+1;
		if (S%xstep==0||S==(s_step + w_step-1)) {ProgressBar ("MCMC sampling ", S, s_step+w_step-1);}
		}
	}

	mixture_no = post_Gn/s_step;
	cout<<endl<<endl<<"MCMC sampling is finished"<<endl;

	VectorXd eigen_alpha = VectorXd::Zero (n_snp);
	eigen_alpha = UtX.transpose()*(bv/s_step)/n_snp;
	WriteCoeff (eigen_alpha, post_Ebeta/s_step);
	pheno_mean = post_Ealpha(0)/s_step;

/////////
	XEbeta.setZero(); //set Gbeta to zero first
		for (size_t i=0; i<n_snp; i++)   {
			x_col   =  UtX.col(i);
			XEbeta += x_col*post_Ebeta(i)/s_step;
		}
	WEalpha.setZero();
		for (size_t j=0;j<n_j; j++)   {
			WEalpha += UtW.col(j)*post_Ealpha(j)/s_step;
		}

	y_res.setZero();
	y_res = Uty - XEbeta - WEalpha;

	double llike_hat = logLike(D,y_res,post_sigma2b/s_step,post_sigma2e/s_step); 
	VectorXd llike2  = llike.tail(s_step).array();
	          llike2 = llike2.array() - llike2.sum()/s_step;
	pD1   =  2*(llike_hat-post_llike/s_step);
	pD2   =  2*(llike2.dot(llike2))/(s_step-1);
	if (pD1 < 0)   {pD1 = 1;}
	DIC1  = -2* llike_hat +  2*pD1;
	DIC2  = -2* llike_hat +  2*pD2;
	BIC1  = -2* llike_hat + log(n_idv)*pD1;
	BIC2  = -2* llike_hat + log(n_idv)*pD2;

	//cout<<" K     = " << n_k<<endl;
	//cout<<" pD1   = " << pD1 <<endl;
	//cout<<" pD2   = " << pD2<<endl;
	//cout<<" DIC1  = " << DIC1<<endl;
	//cout<<" DIC2  = " << DIC2<<endl;
	//cout<<" BIC1  = " << BIC1<<endl;
	//cout<<" BIC2  = " << BIC2<<endl;

	//ofstream outfile1;
	//outfile1.open ("llike.txt", ios::out | ios::binary);
	//for (size_t k = 0; k<s_step; k++)  {
		//outfile1 << setprecision(7) <<llike2(k)<<endl;
		//}
	//outfile1 << endl;


	//cout<<endl<< "Computaion Time for LDR.Gibbs = ";
	//cout<<(clock()-time_begin)/(double(CLOCKS_PER_SEC)*60.0)<<" min"<<endl<<endl;
/*
	//output mik and pik
	ofstream outfile1;
	ofstream outfile2;
	ofstream outfile3;
	ofstream outfile4;
	outfile1.open ("mik.txt", ios::out | ios::binary);
	outfile2.open ("pik.txt", ios::out | ios::binary);
	outfile3.open ("sik.txt", ios::out | ios::binary);
	//outfile4.open ("sigmak.txt", ios::out | ios::binary);
	
	for (size_t i = 0; i<n_snp; i++)  {
		for (size_t k = 0; k<n_k; k++)  {
			outfile1 << setprecision(4) <<  mik_beta(i,k)<<"\t";
			outfile2 << setprecision(4) <<  pik_beta(i,k)<<"\t";
			outfile3 << setprecision(4) << sik2_beta(i,k)<<"\t";
			//outfile4 << setprecision(4) << b_k(k)/a_k(k) <<"\t";
			}
		outfile1 << endl;
		outfile2 << endl;
		outfile3 << endl;
		//outfile4 << endl;
	}
	outfile1.close();
	outfile2.close();
	outfile3.close();
	//outfile4.close();
	//output mik and pik
*/

	////////////
	xtx.resize(0);xty.resize(0);wtw.resize(0);wty.resize(0);
	Utu.resize(0);Ue.resize(0);Ub.resize(0);bv0.resize(0);
	bv.resize(0);V.resize(0);M.resize(0);Ealpha.resize(0);
	post_Ealpha.resize(0);m_alpha.resize(0);s2_alpha.resize(0);
	WEalpha.resize(0);XEbeta.resize(0);Ebeta.resize(0);
	post_Ebeta.resize(0);Ebeta2k.resize(0,0);
	B1.resize(0,0);mik_beta.resize(0,0);sik2_beta.resize(0,0);
	pik_beta.resize(0,0);beta_beta.resize(0,0);gamma_beta.resize(0,0);
	sigma2k.resize(0);vk.resize(0);Elogsigmak.resize(0);
	Elogvk.resize(0);sumElogvl.resize(0);pikexp1.resize(0);
	pikexp2.resize(0);index.resize(0);a_k.resize(0);
	b_k.resize(0);kappa_k.resize(0);lambda_k.resize(0);
	y_res.resize(0);bk.resize(0);
	eigen_alpha.resize(0);
}

void LDR::GibbScreen (const MatrixXd UtX, const VectorXd Uty, const MatrixXd UtW,const VectorXd D, const VectorXd Wbeta, const VectorXd se_Wbeta, const VectorXd beta, const VectorXi snp_no, double lambda) 
{
	
	setprecision(5); 
	clock_t time_begin  = clock();
	size_t T_snp = UtX.cols();
	size_t n_snp = snp_no.size();
	size_t n_idv = Uty.size();
	size_t n_j   = UtW.cols();
	VectorXd x_col(n_idv);
	VectorXd xtx(n_snp); 
	VectorXd xty(n_snp); 
	VectorXd wtw(n_j);		 
	VectorXd wty(n_j);
	VectorXd Utu = VectorXd::Zero(n_idv);
	VectorXd Ue(n_idv);
	VectorXd Ub(n_idv);
	VectorXd bv0 = VectorXd::Zero(n_idv); //for computation of alpha
	VectorXd bv  = VectorXd::Zero(n_idv); //for computation of alpha
	VectorXd V = VectorXd::Zero (n_idv);//compute Utu 
	VectorXd M = VectorXd::Zero (n_idv);//compute Utu

	VectorXd Ealpha      = Wbeta;//intercept
	VectorXd post_Ealpha = VectorXd::Zero(n_j);//save intercept
	VectorXd m_alpha     = Wbeta;
	VectorXd s2_alpha    = se_Wbeta.array()*se_Wbeta.array();

	VectorXd WEalpha(n_idv); 
	double   sigma2e=0, lambdax;
	VectorXd XEbeta(n_idv); //G*beta
	VectorXd Ebeta(n_snp);
	VectorXd post_Ebeta(n_snp);//save beta
	MatrixXd Ebeta2k;
	MatrixXd B1;
	MatrixXd mik_beta = MatrixXd::Zero(n_snp, n_k);

	for (size_t i=0;i<n_snp; i++)  {
		//Eigen_matrix_get_row (UtX, i, x_col);
		x_col  = UtX.col(snp_no(i));
		xtx(i) = x_col.dot(x_col);
		xty(i) = x_col.dot(Uty);
	}

	for (size_t j=0;j<n_j; j++)  {
		wtw(j) = UtW.col(j).dot(UtW.col(j));
		wty(j) = UtW.col(j).dot(Uty);
	}

	//// initial values for SNP and each normal
	//// component has the same initial values
	for (size_t i=0;i<n_snp; i++)  {
		for (size_t k=1;k<n_k; k++)  {
			mik_beta(i,k) = beta(snp_no(i));}}

	mik_beta.col(0).fill(0);
	MatrixXd sik2_beta  =(MatrixXd::Random(n_snp, n_k)).array().abs();
	MatrixXd pik_beta   = MatrixXd::Zero  (n_snp, n_k);
	pik_beta = pik_beta.array() + (double)1/n_k;

	MatrixXd beta_beta  = MatrixXd::Zero(n_snp, n_k);
	MatrixXd gamma_beta = MatrixXd::Zero(n_snp, n_k);
	MatrixXd post_gamma_beta = gamma_beta;
	gamma_beta = gamma_beta.array() + (double)1/n_k;

	VectorXd sigma2k (n_k);
	VectorXd vk (n_k);
	VectorXd Elogsigmak (n_k);
	VectorXd Elogvk (n_k);
	VectorXd sumElogvl (n_k);

	VectorXd pikexp1(n_k);
	VectorXd pikexp2(n_k);
	VectorXd index;
	VectorXd a_k = (VectorXd::Random(n_k)).array().abs();
	VectorXd b_k = (VectorXd::Random(n_k)).array().abs();
	VectorXd  kappa_k = pik_beta.colwise().sum();
	VectorXd lambda_k = (VectorXd::Random(n_k)).array().abs();

	XEbeta.setZero();//set to zero
	for (size_t i=0; i<n_snp; i++)   {
		//Eigen_matrix_get_row (UtX, i, x_col);
		x_col   =  UtX.col(snp_no(i));
		XEbeta += x_col*beta(snp_no(i));
		}

	WEalpha.setZero(); //set to zero
	for (size_t j=0;j<n_j; j++)   {
		WEalpha += UtW.col(j)*m_alpha(j);
		}

	VectorXd y_res  = Uty - WEalpha - XEbeta;
	double sigma2e0 = ((y_res).dot(y_res))/(n_idv-1);
	cout << endl << endl;
	double a_e  = (2*n_idv + n_snp)/2;
	double b_e  = (a_e-1)*sigma2e0;

	double ak   = 21;
	if (lambda < 0.01)     {lambda =    0.01;}
	else if (lambda > 100) {lambda =     100;}
	else                   {lambda =  lambda;}
	//double bk0  = lambda*(ak-1)/n_snp;
	double bk0  = lambda*(ak-1)/UtX.cols();

	VectorXd bk(n_k);
	bk(0) = bk0;
	for (size_t i=1;i<n_k; i++)  {
		bk(i) = bk(i-1)*2.7*sqrt(pow(i,i));
		}

	for (size_t i=0;i<n_snp; i++)  {
		for (size_t k=0;k<n_k; k++)  {
			sik2_beta(i,k) = bk(k)/(ak-1)*sigma2e0;
		}
		}

	double ae   = 0.1;
	double be   = 0.1;
	double a0   = 400;
	double b0   = 40;
	double a_lambda = a0 + n_k;
	double b_lambda = b0;
	double Elogsigmae, tx_ywx_res, xtxabk, A, post_Gn=0;
	double sigma2b = 0.001;

	double post_llike=0, post_sigma2e=0, post_sigma2b=0;

	VectorXd sigma2bX = (VectorXd::Random(w_step+s_step+1)).array().abs();
	sigma2bX(0)=sigma2b;

	VectorXd h = (VectorXd::Random(w_step+s_step+1)).array().abs();
	VectorXd H0(n_idv);
	VectorXd H(n_idv);

	VectorXd    llike = VectorXd::Zero(w_step+s_step);


	////initial values for a_k, b_k and sigma2k
	Ebeta2k = (mik_beta.cwiseProduct(mik_beta)) + sik2_beta;
	for (size_t k=0; k<n_k; k++)  {
		a_k(k) = (pik_beta.col(k).sum())/2 + ak;
		b_k(k) = (Ebeta2k.col(k).sum())*(a_e/b_e)/2 + bk(k);
		sigma2k(k) = b_k(k)/(a_k(k)-1);
		//cout<<sigma2k(k)<<endl;
		}

	A = (y_res).dot(y_res);
	B1=(1/sigma2k.array()).rowwise().replicate(n_snp);
	B1.col(0).fill(0);

	Ebeta2k = Ebeta2k.cwiseProduct(B1.transpose());
	double B = Ebeta2k.sum();
	a_e = n_idv + T_snp*0.01 + ae;
	b_e = (A + B + 2*be)/2;
	////random seed
	gsl_rng * rs = gsl_rng_alloc (gsl_rng_mt19937);

	if (n_k == 1)   {
		WritelmmBeta(beta);
		pheno_mean = m_alpha(0);
		}

	// //  begin MCMC sampling
	else {
	for (size_t S=0; S<(w_step + s_step); S++)   { 

		sigma2e    = 1/gsl_ran_gamma(rs,a_e,1/b_e);
		Elogsigmae = log(sqrt(sigma2e));

		// save Ebeta for the mixture normal component
		Ebeta = (beta_beta.cwiseProduct(gamma_beta)).rowwise().sum();
		
		if (S > (w_step-1))   {post_Ebeta += Ebeta;}

		//sample sigma2k and compute related quantities
		for (size_t k=0; k<n_k; k++)   {
			if (k==0) {
				sigma2k(k)    = 0; 
				Elogsigmak(k) = 0;
				}

			else {
				sigma2k(k)    = 1/gsl_ran_gamma(rs,a_k(k),1/b_k(k));
				Elogsigmak(k) = log(sqrt(sigma2k(k)));
				}

				if (k == (n_k-1)) {vk(k) = 0; Elogvk(k) = 0;}
				else {
					vk(k)     = gsl_ran_beta(rs,kappa_k(k),lambda_k(k));
					Elogvk(k) = log(vk(k));
				}
				sumElogvl(k)  = sum_Elogvl2(vk,k);
			}

		//////////////  sampling the mixture snp effects 
		XEbeta.setZero(); //set Gbeta to zero first
		for (size_t i=0; i<n_snp; i++)   {
			x_col   =  UtX.col(snp_no(i));
			XEbeta += x_col*Ebeta(i);
		}

		y_res.setZero();
		y_res = Uty - WEalpha;

		H0 = (sigma2b*D).array()+1;
		H0 = 1/H0.array();

		//////////////////////////////////////
		for (size_t i=0; i<n_snp; i++)   {

			x_col      = UtX.col(snp_no(i));
			XEbeta    -= x_col*Ebeta(i);

			H = x_col.array()*H0.array();
			tx_ywx_res = H.dot(y_res - XEbeta);

			for (size_t k=0; k<n_k; k++)   {
				
				if (k == 0)   {
					mik_beta(i,k)  = 0; 
					sik2_beta(i,k) = 0; 
					pikexp1(k)     = 0;
					}

				else {
					xtxabk         = H.dot(x_col) + 1/sigma2k(k);
					mik_beta(i,k)  = tx_ywx_res/xtxabk;
					sik2_beta(i,k) = sigma2e/xtxabk;
					pikexp1(k)     = mik_beta(i,k)*mik_beta(i,k)/(2*sik2_beta(i,k)) + log(sqrt(sik2_beta(i,k))) - Elogsigmae - Elogsigmak(k);
				}

				beta_beta(i,k) = gsl_ran_gaussian(rs,sqrt(sik2_beta(i,k))) + mik_beta(i,k);
				pikexp2(k)     = Elogvk(k) + sumElogvl(k);
				}

			index = pikexp1 + pikexp2;
			index = (index.array() - index.maxCoeff()).exp();
			pik_beta.row(i) = index/index.sum();

		// multinomial sampling
		double mult_prob[n_k];
		unsigned int mult_no[n_k];
		for (size_t k=0; k<n_k; k++)   {
			mult_prob[k] = pik_beta(i,k);
			}
		gsl_ran_multinomial(rs, n_k, 1, mult_prob, mult_no);
		for (size_t k=0; k<n_k; k++)   {
			gamma_beta(i,k) = mult_no[k];
		}
		// multinomial sampling

		Ebeta(i) = ((beta_beta.row(i)).dot(gamma_beta.row(i)));
		XEbeta  += x_col*Ebeta(i);
		}
		//////////////////////////////////////

		if (S > (w_step-1))   {
			post_gamma_beta += gamma_beta;
			}
		//////////////  sampling the mixture snp effects 

		WEalpha.setZero();
		for (size_t j=0;j<n_j; j++)   {
			WEalpha += UtW.col(j)*Ealpha(j);
		}
		
		y_res.setZero();
		y_res = Uty - XEbeta;
		for (size_t j=0;j<n_j; j++)   {
			WEalpha    -= UtW.col(j)*Ealpha(j);
					  H = UtW.col(j).array()*H0.array();
			m_alpha(j)  = (H.dot(y_res - WEalpha))/(H.dot(UtW.col(j)));
			s2_alpha(j) = sigma2e/(H.dot(UtW.col(j)));
			Ealpha(j)   = m_alpha(j) + gsl_ran_gaussian(rs,sqrt(s2_alpha(j)));
			WEalpha    += UtW.col(j)*Ealpha(j);
			}

		if (S>(w_step-1))  {post_Ealpha += Ealpha;}

		a_lambda = a0 + n_k;
		b_lambda = b0 - sum_Elogvl2(vk,n_k-1);
		lambdax  = gsl_ran_gamma(rs,a_lambda,1/b_lambda);

		for (size_t k=0; k<n_k-1; k++)   {
			kappa_k(k) = gamma_beta.col(k).sum()+(T_snp-n_snp)*(k==0)+1;
			lambda_k(k) = sum_lambda_k(gamma_beta,k,n_k) + lambdax;
			}

		Ebeta2k = beta_beta.cwiseProduct(beta_beta).cwiseProduct(gamma_beta);
		for (size_t k=0; k<n_k; k++)   {
			a_k(k) = (gamma_beta.col(k).sum()+(T_snp-n_snp)*(k==0))/2 + ak;
			b_k(k) = (Ebeta2k.col(k).sum())/(2*sigma2e) + bk(k);
			}

		y_res.setZero();
		y_res = Uty - XEbeta - WEalpha;

		bv0 = (y_res*sigma2b).array()*H0.array();
		if (S>(w_step-1))   {bv += bv0;}

		H = y_res.array()*H0.array();

		A  = H.dot(y_res);
		B1 = (1/sigma2k.tail(n_k-1).array()).rowwise().replicate(n_snp);
		Ebeta2k   = Ebeta2k.rightCols(n_k-1).cwiseProduct(B1.transpose());
		double B  = Ebeta2k.sum();
		double Gn = gamma_beta.rightCols(n_k-1).sum();
		
		a_e = n_idv/2 + Gn/2 + ae;
		b_e = (A + B + 2*be)/2;

		//double h_new  = gsl_rng_uniform(rs);
		double h_new  = gsl_ran_beta(rs,2,8);
		double ratio1 = log_h(D,y_res,h_new,sigma2e,ae,be)-log(gsl_ran_beta_pdf(h_new,2,8));
		double ratio2 = log_h(D,y_res,h(S), sigma2e,ae,be)-log(gsl_ran_beta_pdf(h(S),2,8));
		double ratio  = exp(ratio1-ratio2);
		if (ratio>1) {ratio = 1;}
		else         {ratio = ratio;}
		double u = gsl_rng_uniform(rs);
		if (u<ratio) {h(S+1) = h_new;}
		else         {h(S+1) = h(S);}

		sigma2b = h(S+1)/(1-h(S+1));

		if (S > (w_step-1)) {post_Gn += Gn;}
		////////////////////////////////////
		double llike0=logLike(D,y_res,sigma2b,sigma2e); 
		     llike(S)=llike0; 

		if (S > (w_step-1))   {	
			post_llike   += llike0;
			post_sigma2e += sigma2e;
			post_sigma2b += sigma2b;
			}

		size_t xstep = (s_step + w_step)*0.2+1;
		if (S%xstep==0||S==(s_step + w_step-1)) {ProgressBar ("screen MCMC sampling ", S, s_step+w_step-1);}
		}
	}
	
	mixture_no = snp_no.size();
	cout<<endl<<endl<<"MCMC sampling is finished and now compute the polygenic effects ..."<<endl;

	VectorXd eigen_alpha = VectorXd::Zero (T_snp);
	VectorXd Beta        = VectorXd::Zero (T_snp);
	eigen_alpha          = UtX.transpose()*(bv/s_step)/T_snp;

	int mcount=0;
	for (size_t i=0; i<T_snp; i++)  {
		if (mcount<snp_no.size())  {
			if (i==(size_t)snp_no(mcount))  {
				Beta(i) = post_Ebeta(mcount)/s_step;
				mcount++;
				}}
		else {Beta(i) = Beta(i);}
	}

	WriteCoeff (eigen_alpha, Beta);
	pheno_mean = post_Ealpha(0)/s_step;

/////////
	XEbeta.setZero(); //set Gbeta to zero first
		for (size_t i=0; i<n_snp; i++)   {
			x_col   =  UtX.col(i);
			XEbeta += x_col*post_Ebeta(i)/s_step;
		}
	WEalpha.setZero();
		for (size_t j=0;j<n_j; j++)   {
			WEalpha += UtW.col(j)*post_Ealpha(j)/s_step;
		}

	y_res.setZero();
	y_res = Uty - XEbeta - WEalpha;

	double llike_hat = logLike(D,y_res,post_sigma2b/s_step,post_sigma2e/s_step); 
	VectorXd llike2  = llike.tail(s_step).array();
	          llike2 = llike2.array() - llike2.sum()/s_step;
	pD1   =  2*(llike_hat-post_llike/s_step);
	pD2   =  2*(llike2.dot(llike2))/(s_step-1);
	if (pD1 < 0)   {pD1 = 1;}
	DIC1  = -2* llike_hat +  2*pD1;
	DIC2  = -2* llike_hat +  2*pD2;
	BIC1  = -2* llike_hat + log(n_idv)*pD1;
	BIC2  = -2* llike_hat + log(n_idv)*pD2;

	cout<<" K     = " << n_k<<endl;
	cout<<" pD1   = " << pD1 <<endl;
	cout<<" pD2   = " << pD2<<endl;
	cout<<" DIC1  = " << DIC1<<endl;
	cout<<" DIC2  = " << DIC2<<endl;
	cout<<" BIC1  = " << BIC1<<endl;
	cout<<" BIC2  = " << BIC2<<endl;

	cout<<endl<< "Computaion Time for LDR.Gibbs = ";
	cout<<(clock()-time_begin)/(double(CLOCKS_PER_SEC)*60.0)<<" min"<<endl<<endl;

	////////////
	xtx.resize(0);xty.resize(0);wtw.resize(0);wty.resize(0);Utu.resize(0);
	Ue.resize(0);Ub.resize(0);bv0.resize(0);bv.resize(0);V.resize(0);
	M.resize(0);Ealpha.resize(0);post_Ealpha.resize(0);m_alpha.resize(0);
	s2_alpha.resize(0);WEalpha.resize(0);XEbeta.resize(0);Ebeta.resize(0);
	post_Ebeta.resize(0);Ebeta2k.resize(0,0);B1.resize(0,0);mik_beta.resize(0,0);
	sik2_beta.resize(0,0);pik_beta.resize(0,0);beta_beta.resize(0,0);
	gamma_beta.resize(0,0);sigma2k.resize(0);vk.resize(0);Elogsigmak.resize(0);
	Elogvk.resize(0);sumElogvl.resize(0);pikexp1.resize(0);pikexp2.resize(0);
	index.resize(0);a_k.resize(0);b_k.resize(0);kappa_k.resize(0);
	lambda_k.resize(0);y_res.resize(0);bk.resize(0);Beta.resize(0);
}



void LDR::GibbScreenX (const MatrixXd UtX, const VectorXd Uty, const MatrixXd UtW,const VectorXd D, const VectorXd Wbeta, const VectorXd se_Wbeta, const VectorXd beta, const VectorXi snp_no, double lambda) 
{
	
	setprecision(5); 
	clock_t time_begin  = clock();
	//size_t n_snp = UtX.cols();
	size_t n_snp = snp_no.size();
	size_t n_idv = Uty.size();
	size_t n_j   = UtW.cols();
	VectorXd x_col(n_idv);
	VectorXd xtx(n_snp); 
	VectorXd xty(n_snp); 
	VectorXd wtw(n_j);		 
	VectorXd wty(n_j);
	VectorXd Utu = VectorXd::Zero(n_idv);
	VectorXd Ue(n_idv);
	VectorXd Ub(n_idv);
	VectorXd bv0 = VectorXd::Zero(n_idv); //for computation of alpha
	VectorXd bv  = VectorXd::Zero(n_idv); //for computation of alpha
	VectorXd V = VectorXd::Zero (n_idv);//compute Utu 
	VectorXd M = VectorXd::Zero (n_idv);//compute Utu
	VectorXd snp_label(n_snp);

	VectorXd Ealpha      = Wbeta;//intercept
	VectorXd post_Ealpha = VectorXd::Zero(n_j);//save intercept
	VectorXd m_alpha     = Wbeta;
	VectorXd s2_alpha    = se_Wbeta.array()*se_Wbeta.array();

	VectorXd WEalpha(n_idv); 
	double   sigma2e=0, lambdax;
	VectorXd XEbeta(n_idv); //G*beta
	VectorXd Ebeta(n_snp);
	VectorXd post_Ebeta(n_snp);//save beta
	MatrixXd Ebeta2k;
	MatrixXd B1;
	MatrixXd mik_beta = MatrixXd::Zero(n_snp, n_k);

	for (size_t i=0;i<n_snp; i++)  {
		//Eigen_matrix_get_row (UtX, i, x_col);
		x_col  = UtX.col(snp_no(i));
		xtx(i) = x_col.dot(x_col);
		xty(i) = x_col.dot(Uty);
	}

	for (size_t j=0;j<n_j; j++)  {
		wtw(j) = UtW.col(j).dot(UtW.col(j));
		wty(j) = UtW.col(j).dot(Uty);
	}

	//// initial values for SNP and each normal
	//// component has the same initial values
	for (size_t i=0;i<n_snp; i++)  {
		for (size_t k=1;k<n_k; k++)  {
			mik_beta(i,k) = beta(snp_no(i));
		}
	}

	mik_beta.col(0).fill(0);
	MatrixXd sik2_beta  =(MatrixXd::Random(n_snp, n_k)).array().abs();
	MatrixXd pik_beta   = MatrixXd::Zero  (n_snp, n_k);
	pik_beta = pik_beta.array() + (double)1/n_k;

	//MatrixXd beta_beta  = mik_beta;
	MatrixXd beta_beta  = MatrixXd::Zero(n_snp, n_k);
	MatrixXd gamma_beta = MatrixXd::Zero(n_snp, n_k);
	MatrixXd post_gamma_beta = gamma_beta;
	MatrixXd post_pik_beta   = gamma_beta;
	gamma_beta = gamma_beta.array() + (double)1/n_k;

	VectorXd sigma2k (n_k);
	VectorXd vk (n_k);
	VectorXd Elogsigmak (n_k);
	VectorXd Elogvk (n_k);
	VectorXd sumElogvl (n_k);

	VectorXd pikexp1(n_k);
	VectorXd pikexp2(n_k);
	VectorXd index;
	VectorXd a_k = (VectorXd::Random(n_k)).array().abs();
	VectorXd b_k = (VectorXd::Random(n_k)).array().abs();
	VectorXd  kappa_k = pik_beta.colwise().sum();
	VectorXd lambda_k = (VectorXd::Random(n_k)).array().abs();

	XEbeta.setZero();//set to zero
	for (size_t i=0; i<n_snp; i++)   {
		//Eigen_matrix_get_row (UtX, i, x_col);
		x_col   =  UtX.col(snp_no(i));
		XEbeta += x_col*beta(snp_no(i));
		}

	WEalpha.setZero(); //set to zero
	for (size_t j=0;j<n_j; j++)   {
		WEalpha += UtW.col(j)*m_alpha(j);
		}

	VectorXd y_res  = Uty - WEalpha - XEbeta;
	double sigma2e0 = ((y_res).dot(y_res))/(n_idv-1);
	//cout << endl << "sigma2e0 = " << sigma2e0 << endl;
	cout << endl << endl;
	double a_e  = (2*n_idv + n_snp)/2;
	double b_e  = (a_e-1)*sigma2e0;

	double ak   = 11;
	if (lambda < 0.01)     {lambda =    0.01;}
	else if (lambda > 100) {lambda =     100;}
	else                   {lambda =  lambda;}
	double bk0  = lambda*(ak-1)/n_snp;

	VectorXd bk(n_k);
	bk(0) = bk0;
	for (size_t i=1;i<n_k; i++)  {
		bk(i) = bk(i-1)*1.7*sqrt(pow(i,i));
		}

	for (size_t i=0;i<n_snp; i++)  {
		for (size_t k=0;k<n_k; k++)  {
			sik2_beta(i,k) = bk(k)/(ak-1)*sigma2e0;
			}
		}

	double ae   = 0.1;
	double be   = 0.1;
	double a0   = 400;
	double b0   = 40;
	double a_lambda = a0 + n_k;
	double b_lambda = b0;
	double Elogsigmae, tx_ywx_res, xtxabk, A;
	double sigma2b = 0.001;

	////initial values for a_k, b_k and sigma2k
	Ebeta2k = (mik_beta.cwiseProduct(mik_beta)) + sik2_beta;
	for (size_t k=0; k<n_k; k++)  {
		a_k(k) = (pik_beta.col(k).sum())/2 + ak;
		b_k(k) = (Ebeta2k.col(k).sum())*(a_e/b_e)/2 + bk(k);
		sigma2k(k) = b_k(k)/(a_k(k)-1);
		//cout<<sigma2k(k)<<endl;
		}

	A = (y_res).dot(y_res);
	B1=(1/sigma2k.array()).rowwise().replicate(n_snp);
	//B1.col(0).fill(0);

	Ebeta2k = Ebeta2k.cwiseProduct(B1.transpose());
	double B = Ebeta2k.sum();
	a_e = n_idv + n_snp*0.1 + ae;
	b_e = (A + B + 2*be)/2;
	////random seed
	gsl_rng * rs = gsl_rng_alloc (gsl_rng_mt19937);

	//WritelmmBeta(beta);
	if (n_k == 1)   {
		WritegibbsBeta (beta);
		pheno_mean = m_alpha(0);
		}

	// //  begin MCMC sampling
	else {
	for (size_t S=0; S<(w_step + s_step); S++)   { 

		sigma2e    = 1/gsl_ran_gamma(rs,a_e,1/b_e);
		Elogsigmae = log(sqrt(sigma2e));

		// save Ebeta for the mixture normal component
		Ebeta = (beta_beta.cwiseProduct(gamma_beta)).rowwise().sum();		
		if (S > (w_step-1))   {
			post_Ebeta += Ebeta;
			}

		//sample sigma2k and compute related quantities
		for (size_t k=0; k<n_k; k++)   {

			//if (k==0) {
				//sigma2k(k)    = 0; 
				//Elogsigmak(k) = 0;
				//}

			//else {
				sigma2k(k)    = 1/gsl_ran_gamma(rs,a_k(k),1/b_k(k));
				Elogsigmak(k) = log(sqrt(sigma2k(k)));
				//}

				if (k == (n_k-1)) {vk(k) = 0; Elogvk(k) = 0;}
				else {
					vk(k)     = gsl_ran_beta(rs,kappa_k(k),lambda_k(k));
					Elogvk(k) = log(vk(k));
				}
				sumElogvl(k)  = sum_Elogvl2(vk,k);
			}

		//////////////  sampling the mixture snp effects 
		XEbeta.setZero(); //set Gbeta to zero first
			for (size_t i=0; i<n_snp; i++)   {
				//Eigen_matrix_get_row (UtX, i, x_col);
				x_col   =  UtX.col(snp_no(i));
				XEbeta += x_col*Ebeta(i);
		}

		y_res.setZero();
		y_res = Uty - WEalpha - Utu;
	
		//////////////////////////////////////
		for (size_t i=0; i<n_snp; i++)   {

			//Eigen_matrix_get_row (UtX, i, x_col);
			x_col      = UtX.col(snp_no(i));
			XEbeta    -= x_col*Ebeta(i);
			tx_ywx_res = x_col.dot(y_res - XEbeta);
			for (size_t k=0; k<n_k; k++)   {				
				//if (k == 0)   {
					//mik_beta(i,k)  = 0; 
					//sik2_beta(i,k) = 0; 
					//pikexp1(k)     = 0;
					//}
				//else {
					xtxabk         = xtx(i) + 1/sigma2k(k);
					mik_beta(i,k)  = tx_ywx_res/xtxabk;
					sik2_beta(i,k) = sigma2e/xtxabk;
					pikexp1(k)     = mik_beta(i,k)*mik_beta(i,k)/(2*sik2_beta(i,k)) + log(sqrt(sik2_beta(i,k))) - Elogsigmae - Elogsigmak(k);
				//}
				beta_beta(i,k) = gsl_ran_gaussian(rs,sqrt(sik2_beta(i,k))) + mik_beta(i,k);
				pikexp2(k)     = Elogvk(k) + sumElogvl(k);
				}
			index = pikexp1 + pikexp2;
			index = (index.array() - index.maxCoeff()).exp();
			pik_beta.row(i) = index/index.sum();

		// multinomial sampling
		double mult_prob[n_k];
		unsigned int mult_no[n_k];
		for (size_t k=0; k<n_k; k++)   {
			mult_prob[k] = pik_beta(i,k);
			}
		gsl_ran_multinomial(rs, n_k, 1, mult_prob, mult_no);
		for (size_t k=0; k<n_k; k++)   {
			gamma_beta(i,k) = mult_no[k];
		}
		// multinomial sampling

		//non_zero+= gamma_beta.row(i).tail(n_k-1).sum();
		Ebeta(i) = ((beta_beta.row(i)).dot(gamma_beta.row(i)));
		XEbeta  += x_col*Ebeta(i);
		//cout<<pik_beta.row(i)<<endl;
		}
		//////////////////////////////////////

		if (S > (w_step-1))   {
			post_gamma_beta += gamma_beta;
			post_pik_beta   += pik_beta;
			}
		//////////////  sampling the mixture snp effects 

		WEalpha.setZero();
		for (size_t j=0;j<n_j; j++)   {
			WEalpha += UtW.col(j)*Ealpha(j);
		}
		
		y_res.setZero();
		y_res = Uty - XEbeta - Utu;
		for (size_t j=0;j<n_j; j++)   {
			WEalpha    -= UtW.col(j)*Ealpha(j);
			m_alpha(j)  = 1/wtw(j)*(UtW.col(j).dot(y_res - WEalpha));
			s2_alpha(j) = sigma2e/wtw(j);
			Ealpha(j)   = m_alpha(j) + gsl_ran_gaussian(rs,sqrt(s2_alpha(j)));
			WEalpha    += UtW.col(j)*Ealpha(j);
			}

		if (S>(w_step-1))  {post_Ealpha += Ealpha;}

		a_lambda = a0 + n_k;
		b_lambda = b0 - sum_Elogvl2(vk,n_k-1);
		lambdax  = gsl_ran_gamma(rs,a_lambda,1/b_lambda);

		for (size_t k=0; k<n_k-1; k++)   {
			kappa_k(k)  = gamma_beta.col(k).sum() + 1;
			lambda_k(k) = sum_lambda_k(gamma_beta,k,n_k) + lambdax;
			}

		Ebeta2k = beta_beta.cwiseProduct(beta_beta).cwiseProduct(gamma_beta);
		for (size_t k=0; k<n_k; k++)   {
			a_k(k) = (gamma_beta.col(k).sum())/2 + ak;
			b_k(k) = (Ebeta2k.col(k).sum())/(2*sigma2e) + bk(k);
			}

		y_res.setZero();
		y_res = Uty - XEbeta - WEalpha;
		for (size_t i=0; i<n_idv; i++) {
			V(i)  = sigma2b*D(i)/(sigma2b*D(i) + 1);
			M(i)  = y_res(i)*V(i);
			Utu(i)= M(i) + gsl_ran_gaussian(rs,sqrt(V(i)*sigma2e));
			if (D(i)  == 0) {
				Ue(i) = 0; 
				Ub(i) = 0;
				}
			else {
				Ue(i) = Utu(i)/(sigma2b*D(i)); //for sigma2e
				Ub(i) = Utu(i)/(sigma2e*D(i)); //for sigma2b
				}
			bv0(i) = y_res(i)*sigma2b/(sigma2b*D(i) + 1); 
			}

		if (S>(w_step-1)) {bv += bv0;}

		A  = (y_res-Utu).dot(y_res-Utu);
		//B1 = (1/sigma2k.tail(n_k-1).array()).rowwise().replicate(n_snp);
		B1=(1/sigma2k.array()).rowwise().replicate(n_snp);
		Ebeta2k = Ebeta2k.cwiseProduct(B1.transpose());
		//Ebeta2k   = Ebeta2k.rightCols(n_k-1).cwiseProduct(B1.transpose());
		double B  = Ebeta2k.sum();
		//double Gn = gamma_beta.rightCols(n_k-1).sum();
		double Gn = gamma_beta.sum();
		
		a_e = n_idv + Gn/2 + ae;
		b_e = (A + B + Ue.dot(Utu) + 2*be)/2;
		sigma2b = 1/gsl_ran_gamma(rs,n_idv/2+ae,1/(2*Ub.dot(Utu))+be);

		//if (S > (w_step-1)) {post_Gn += Gn;}
		////////////////////////////////////

		//cout<<S+1<<"  "<<setprecision(0)<<Gn<<"  "<<setprecision(4)<<sigma2e<<endl;
		size_t xstep = (s_step + w_step)*0.2;
		if (S%xstep==0||S==(s_step + w_step-1)) {ProgressBar ("MCMC sampling ", S, s_step+w_step-1);}
		//if (S%1==0||S==(s_step + w_step-1)) {ProgressBar ("MCMC sampling ", S, s_step+w_step-1);}
		}
	}
	
	//mixture_no = post_Gn/s_step;
	mixture_no = snp_no.size();
	cout<<endl<<endl<<"MCMC sampling is finished";
	cout<<" and now compute the polygenic effects ..."<<endl;

	size_t T_snp  = UtX.cols();
	VectorXd Alpha = VectorXd::Zero (T_snp);
	VectorXd Beta  = VectorXd::Zero (T_snp);
	Alpha = UtX.transpose()*(bv/s_step)/T_snp;

	int mcount=0;
	for (size_t i=0; i<T_snp; i++)  {
		if (mcount<snp_no.size())  {
			if (i==(size_t)snp_no(mcount))  {
				Beta(i) = post_Ebeta(mcount)/s_step;
				mcount++;
				}}
		else {Beta(i) = Beta(i);}
	}

	WriteCoeff (Alpha, Beta);
	pheno_mean = post_Ealpha(0)/s_step;

	cout << endl;
	cout << "Computaion Time for LDR.Gibbs = ";
	cout << (clock()-time_begin)/(double(CLOCKS_PER_SEC)*60.0)<<" min"<<endl;
	cout << endl;


	//output mik and pik
	ofstream outfile1;
	ofstream outfile2;
	ofstream outfile3;
	ofstream outfile4;
	outfile1.open ("mik.txt", ios::out | ios::binary);
	outfile2.open ("pik.txt", ios::out | ios::binary);
	outfile3.open ("sik.txt", ios::out | ios::binary);
	//outfile4.open ("sigmak.txt", ios::out | ios::binary);
	
	for (size_t i = 0; i<n_snp; i++)  {
		for (size_t k = 0; k<n_k; k++)  {
			outfile1 << setprecision(4) <<  mik_beta(i,k)<<"\t";
			outfile2 << setprecision(4) <<  pik_beta(i,k)<<"\t";
			outfile3 << setprecision(4) << sik2_beta(i,k)<<"\t";
			//outfile4 << setprecision(4) << b_k(k)/a_k(k) <<"\t";
			}
		outfile1 << endl;
		outfile2 << endl;
		outfile3 << endl;
		//outfile4 << endl;
	}
	outfile1.close();
	outfile2.close();
	outfile3.close();
	//outfile4.close();
	//output mik and pik

	////////////
	xtx.resize(0);xty.resize(0);wtw.resize(0);wty.resize(0);Utu.resize(0);
	Ue.resize(0);Ub.resize(0);bv0.resize(0);bv.resize(0);V.resize(0);
	M.resize(0);Ealpha.resize(0);post_Ealpha.resize(0);m_alpha.resize(0);
	s2_alpha.resize(0);WEalpha.resize(0);XEbeta.resize(0);Ebeta.resize(0);
	post_Ebeta.resize(0);Ebeta2k.resize(0,0);B1.resize(0,0);mik_beta.resize(0,0);
	sik2_beta.resize(0,0);pik_beta.resize(0,0);beta_beta.resize(0,0);
	gamma_beta.resize(0,0);sigma2k.resize(0);vk.resize(0);Elogsigmak.resize(0);
	Elogvk.resize(0);sumElogvl.resize(0);pikexp1.resize(0);pikexp2.resize(0);
	index.resize(0);a_k.resize(0);b_k.resize(0);kappa_k.resize(0);
	lambda_k.resize(0);y_res.resize(0);bk.resize(0);Alpha.resize(0);Beta.resize(0);
}


void LDR::VBslmm (const MatrixXd UtX, const VectorXd Uty, const MatrixXd UtW,const VectorXd D, const VectorXd Wbeta, const VectorXd se_Wbeta, const VectorXd beta, double lambda) 
{
	
	setprecision(5); 
	clock_t time_begin  = clock();
	size_t n_snp = UtX.cols();
	size_t n_idv = Uty.size();
	size_t n_j   = UtW.cols();
	VectorXd x_col(n_idv);

	VectorXd xtx(n_snp); 
	VectorXd xty(n_snp); 
	VectorXd wtw(n_j); 
	VectorXd wty(n_j);
	VectorXd g = VectorXd::Zero(n_idv);
	VectorXd Ue(n_idv);
	VectorXd Ub(n_idv);
	VectorXd bv  = VectorXd::Zero(n_idv); //for computation of alpha
	VectorXd V = VectorXd::Zero (n_idv);//compute g 
	VectorXd M = VectorXd::Zero (n_idv);//compute g
	VectorXd snp_label(n_snp);

	VectorXd Ealpha      = Wbeta;//intercept
	VectorXd post_Ealpha = VectorXd::Zero(n_j);//save intercept
	VectorXd m_alpha     = Wbeta;
	VectorXd s2_alpha    = se_Wbeta.array()*se_Wbeta.array();

	VectorXd WEalpha(n_idv); 
	VectorXd XEbeta(n_idv); //G*beta
	VectorXd Ebeta(n_snp);
	VectorXd post_Ebeta(n_snp);//save beta
	VectorXd Ebeta2k(n_snp);

	VectorXd mik_beta = VectorXd::Zero(n_snp);
	mik_beta = beta;

	for (size_t i=0;i<n_snp; i++)  {
		x_col  = UtX.col(i);
		xtx(i) = x_col.dot(x_col);
		xty(i) = x_col.dot(Uty);
	}

	for (size_t j=0;j<n_j; j++)  {
		wtw(j) = UtW.col(j).dot(UtW.col(j));
		wty(j) = UtW.col(j).dot(Uty);
	}

	VectorXd sik2_beta  =(VectorXd::Random(n_snp)).array().abs();
	VectorXd pik_beta   = VectorXd::Zero  (n_snp);
	pik_beta = pik_beta.array() + 0.01;

	VectorXd beta_beta  = VectorXd::Zero(n_snp);
	VectorXd gamma_beta = VectorXd::Zero(n_snp);
	VectorXd post_gamma_beta = gamma_beta;
	VectorXd post_pik_beta   = gamma_beta;
	gamma_beta = gamma_beta.array() + 0.01;

	double Elogsigmak, pikexp1, a_k, b_k;
	double a_pai = n_snp*0.01;
	double b_pai = n_snp-a_pai;

	XEbeta.setZero();//set to zero
	for (size_t i=0; i<n_snp; i++)   {
		x_col   =  UtX.col(i);
		XEbeta += x_col*beta(i);
		}

	WEalpha.setZero(); //set to zero
	for (size_t j=0;j<n_j; j++)   {
		WEalpha += UtW.col(j)*m_alpha(j);
		}

	VectorXd y_res  = Uty - WEalpha - XEbeta;
	double sigma2e0 = ((y_res).dot(y_res))/(n_idv-1);
	cout << endl << "sigma2e0 = " << sigma2e0 << endl;
	cout << endl;
	double a_e  = n_idv + n_snp*0.001/2;
	double b_e  = a_e;

	double ak   = 0.1;
	if (lambda < 0.01)     {lambda =    0.01;}
	else if (lambda > 100) {lambda =     100;}
	else                   {lambda =  lambda;}

	double bk = ak;

	for (size_t i=0;i<n_snp; i++)  {
			sik2_beta(i) = 1/xtx(i);
		}

	double ae   = 0.1;
	double be   = 0.1;
	double a_b  = 0.1;
	double b_b  = 0.1;
	double Elogsigmae, tx_ywx_res, xtxabk, A;

	////initial values for a_k, b_k and sigma2k
	Ebeta2k = mik_beta.array()*mik_beta.array() + sik2_beta.array();
	Ebeta   = mik_beta.array()*pik_beta.array();
	a_k = pik_beta.sum()/2 + ak;
	b_k = Ebeta2k.sum()*(a_e/b_e)/2 + bk;

	size_t int_step = 0;
	size_t max_step = n_idv*sqrt(10);
	double delta = 10;
	VectorXd ELBO = VectorXd::Zero(max_step);
	
	//WritelmmBeta(beta);
	if (n_k == 1)   {
		WritelmmBeta(beta);
		pheno_mean = m_alpha(0);
		}

	/////////////////// variational bayesian/////////////
	else {
	while ((int_step<max_step) && (delta > 1e-5)) {

		Elogsigmae = 0.5*(log(b_e)-gsl_sf_psi(a_e));
		cout<<setprecision(8)<<"Elogsigmae "<<Elogsigmae<<endl;
		//////////////  sampling the mixture snp effects 
		XEbeta.setZero(); //set Gbeta to zero first
		for (size_t i=0; i<n_snp; i++)   {
			x_col   =  UtX.col(i);
			XEbeta += x_col*Ebeta(i);
		}

		//snp_label = permute(n_snp);
		double ab_e = a_e/b_e; // E(sigma2e-2)
		double ba_e = b_e/a_e; // 1/E(sigma2e-2)
		cout<<setprecision(8)<<"ab_e ba_e "<<ab_e<<" "<<ba_e<<endl;

		Elogsigmak     = 0.5*(log(b_k)-gsl_sf_psi(a_k));
		double pikexp2 = -Elogsigmae-Elogsigmak+gsl_sf_psi(a_pai)-gsl_sf_psi(b_pai);
		cout<<setprecision(8)<<"pikexp2 "<<pikexp2<<endl;

		y_res.setZero();
		y_res = Uty - WEalpha - g;
		for (size_t i=0; i<n_snp; i++)   {
			x_col      = UtX.col(i);
			XEbeta    -= x_col*Ebeta(i);
			tx_ywx_res = x_col.dot(y_res - XEbeta);

			xtxabk       = xtx(i) + a_k/b_k;
			mik_beta(i)  = tx_ywx_res/xtxabk;
			sik2_beta(i) = ba_e/xtxabk;

			pikexp1      = pow(mik_beta(i),2)/(2*sik2_beta(i))+0.5*log(sik2_beta(i));
			pik_beta(i)  = exp(pikexp1+pikexp2)/(1+exp(pikexp1+pikexp2));

			Ebeta(i) = mik_beta(i)*pik_beta(i);
			XEbeta  += x_col*Ebeta(i);
			//cout<<setprecision(8)<<"pik mik sik2 "<<pik_beta(i)<<" "<<mik_beta(i)<<" "<<sik2_beta(i)<<endl;
		}
		//////////////  sampling the mixture snp effects 

		WEalpha.setZero();
		for (size_t j=0;j<n_j; j++)   {
			WEalpha += UtW.col(j)*Ealpha(j);
		}
		
		y_res.setZero();
		y_res = Uty - XEbeta - g;
		for (size_t j=0;j<n_j; j++)   {
			WEalpha    -= UtW.col(j)*Ealpha(j);
			Ealpha(j)   = 1/wtw(j)*(UtW.col(j).dot(y_res - WEalpha));
			s2_alpha(j) = ba_e/wtw(j);
			WEalpha    += UtW.col(j)*Ealpha(j);
			cout<<setprecision(8)<<"alpha and s2 "<<Ealpha(j)<<" "<<s2_alpha(j)<<endl;
			}

		a_pai = pik_beta.sum()+1;
		b_pai = n_snp - a_pai+1;

		cout<<setprecision(8)<<"pai "<<a_pai<<" "<<b_pai<<endl;

		Ebeta2k = mik_beta.array()*mik_beta.array() + sik2_beta.array();
		Ebeta2k =  Ebeta2k.array()*pik_beta.array();

		a_k = pik_beta.sum()/2 + ak;
		b_k = Ebeta2k.sum()*ab_e/2 + bk;
		cout<<setprecision(8)<<"a_k  b_k "<<a_k<<" "<<b_k<<endl;

		y_res.setZero();
		y_res = Uty - XEbeta - WEalpha;
		double ab = a_b/b_b;
		for (size_t i=0; i<n_idv; i++)   {
			if (D(i)  == 0)   {
				V(i)  = 0;
				g(i)= 0;
				Ue(i) = 0; 
				Ub(i) = 0;
				}
			else   {
				double abD = ab/D(i);
				V(i)  = ba_e/(abD + 1);
				g(i)= y_res(i)/(abD + 1);
				Ue(i) = (g(i)*g(i)+V(i))*abD;       //for sigma2e
				Ub(i) = (g(i)*g(i)+V(i))*ab_e/D(i); //for sigma2b
			}
			}

		VectorXd VarBeta = Ebeta2k-Ebeta.cwiseProduct(Ebeta);

		A =(y_res-g).dot(y_res-g)+wtw.dot(s2_alpha)+V.sum()+xtx.dot(VarBeta);
		double B  = (a_k/b_k)*Ebeta2k.sum();
		double Gn = pik_beta.sum();
		cout<<setprecision(8)<<"Gn "<<Gn<<endl;

		a_e = n_idv + Gn/2 + ae;
		b_e = (A + B + Ue.sum())/2 + be;

		a_b = n_idv/2 + ae;
		b_b = Ub.sum()/2 + be;

		cout<<setprecision(8)<<"a_e  b_e a_b b_b "<<a_e<<" "<<b_e;
		cout<<" "<<a_b<<" "<<b_b<<endl;

		int_step++;

		ELBO(int_step)= lgamma(a_e)-a_e*log(b_e)  +
						lgamma(a_b)-a_b*log(b_b)  + a_b + 
						lgamma(a_k)-a_e*log(b_k)  + a_k +
						lgamma(a_pai)+lgamma(b_pai)-lgamma(a_pai+b_pai) -
						ELBO3(pik_beta,sik2_beta) +  
						0.5*( s2_alpha.array().log().sum()) + 
						0.5*((V.array()+1e-10).log().sum()) -
						be*(a_k/b_k+ab);

		////////////////////////////////////
		delta = abs((ELBO(int_step)-ELBO(int_step-1))/ELBO(int_step));
		
		if ((int_step+1)%1==0) {cout<<int_step+1<<" "<<setprecision(5)<<delta<<" "<<ELBO(int_step)<<endl;}
		}
	}

	y_res.setZero();
	y_res = Uty - XEbeta - WEalpha;
	for (size_t i=0; i<n_idv; i++)   {
		bv(i) = y_res(i)/(D(i) + a_b/b_b); 
	}

	cout<<endl<<"variational bayes is finished";
	cout<<" and now compute the polygenic effects ..."<<endl;
	VectorXd eigen_alpha = VectorXd::Zero (n_snp);
	eigen_alpha = UtX.transpose()*bv/n_snp;
	//WritevbBeta (Ebeta+eigen_alpha);
	WriteCoeff(eigen_alpha, Ebeta);
	pheno_mean = Ealpha(0);

	cout << endl;
	cout << "delta = "<<scientific<<setprecision(6)<<delta << endl;
	cout << "iteraion count = "<<scientific<<int_step+1<<endl;
	cout << "ELBO = "<<scientific<<setprecision(6)<<ELBO(int_step) << endl;

	cout << endl;
	cout << "Computaion Time for LDR.VB = ";
	cout << setprecision(3)<<(clock()-time_begin)/(double(CLOCKS_PER_SEC)*60.0);
	cout <<" min"<<endl;
	cout << endl;

	////////////
	/*xtx.resize(0);xty.resize(0);wtw.resize(0);wty.resize(0);
	Utu.resize(0);Ue.resize(0);Ub.resize(0);bv.resize(0);
	V.resize(0);M.resize(0);Ealpha.resize(0);post_Ealpha.resize(0);
	m_alpha.resize(0);s2_alpha.resize(0);WEalpha.resize(0);
	XEbeta.resize(0);Ebeta.resize(0);post_Ebeta.resize(0);
	Ebeta2k.resize(0,0);B1.resize(0,0);mik_beta.resize(0,0);
	sik2_beta.resize(0,0);pik_beta.resize(0,0);beta_beta.resize(0,0);
	gamma_beta.resize(0,0);sigma2k.resize(0);vk.resize(0);
	Elogsigmak.resize(0);Elogvk.resize(0);sumElogvl.resize(0);
	pikexp1.resize(0);pikexp2.resize(0);index.resize(0);
	a_k.resize(0);b_k.resize(0);kappa_k.resize(0);
	lambda_k.resize(0);y_res.resize(0);bk.resize(0);
	eigen_alpha.resize(0);*/
}



void LDR::GibbScreen500 (const MatrixXd UtX, const VectorXd Uty, const MatrixXd UtW,const VectorXd D, const VectorXd Wbeta, const VectorXd se_Wbeta, const VectorXd beta, const VectorXi snp_no, double lambda) 
{
	
	setprecision(5); 
	clock_t time_begin  = clock();
	size_t T_snp = UtX.cols();
	//size_t n_snp = snp_no.size();
	size_t n_idv = Uty.size();
	size_t n_j   = UtW.cols();
	VectorXd x_col(n_idv);

	//VectorXd xtx(n_snp); 
	//VectorXd xty(n_snp); 
	VectorXd xtx(T_snp); 
	VectorXd xty(T_snp); 

	VectorXd wtw(n_j);		 
	VectorXd wty(n_j);
	VectorXd Utu = VectorXd::Zero(n_idv);
	VectorXd Ue(n_idv);
	VectorXd Ub(n_idv);
	VectorXd bv0 = VectorXd::Zero(n_idv); //for computation of alpha
	VectorXd bv  = VectorXd::Zero(n_idv); //for computation of alpha
	VectorXd V = VectorXd::Zero (n_idv);//compute Utu 
	VectorXd M = VectorXd::Zero (n_idv);//compute Utu

	VectorXd Ealpha      = Wbeta;//intercept
	VectorXd post_Ealpha = VectorXd::Zero(n_j);//save intercept
	VectorXd m_alpha     = Wbeta;
	VectorXd s2_alpha    = se_Wbeta.array()*se_Wbeta.array();

	VectorXd WEalpha(n_idv); 
	double   sigma2e=0, lambdax;
	VectorXd XEbeta(n_idv); //G*beta

	//VectorXd Ebeta(n_snp);
	//VectorXd post_Ebeta(n_snp);//save beta
	//MatrixXd mik_beta = MatrixXd::Zero(n_snp, n_k);

	VectorXd Ebeta(T_snp);
	VectorXd post_Ebeta(T_snp);//save beta
	MatrixXd mik_beta = MatrixXd::Zero(T_snp, n_k);
		
	MatrixXd Ebeta2k;
	MatrixXd B1;


	//for (size_t i=0;i<n_snp; i++)  {
	for (size_t i=0;i<T_snp; i++)  {
		//Eigen_matrix_get_row (UtX, i, x_col);
		//x_col  = UtX.col(snp_no(i));
		x_col  = UtX.col(i);
		xtx(i) = x_col.dot(x_col);
		xty(i) = x_col.dot(Uty);
	}

	for (size_t j=0;j<n_j; j++)  {
		wtw(j) = UtW.col(j).dot(UtW.col(j));
		wty(j) = UtW.col(j).dot(Uty);
	}

	//// initial values for SNP and each normal
	//// component has the same initial values
	//for (size_t i=0;i<n_snp; i++)  {
	for (size_t i=0;i<T_snp; i++)  {
		for (size_t k=1;k<n_k; k++)  {
			//mik_beta(i,k) = beta(snp_no(i));}}
			mik_beta(i,k) = beta(i);}}

	mik_beta.col(0).fill(0);

	//MatrixXd sik2_beta  =(MatrixXd::Random(n_snp, n_k)).array().abs();
	//MatrixXd pik_beta   = MatrixXd::Zero  (n_snp, n_k);
	//MatrixXd beta_beta  = MatrixXd::Zero(n_snp, n_k);
	//MatrixXd gamma_beta = MatrixXd::Zero(n_snp, n_k);

	MatrixXd sik2_beta  =(MatrixXd::Random(T_snp, n_k)).array().abs();
	MatrixXd pik_beta   = MatrixXd::Zero  (T_snp, n_k);
	MatrixXd beta_beta  = MatrixXd::Zero(T_snp, n_k);
	MatrixXd gamma_beta = MatrixXd::Zero(T_snp, n_k);

	pik_beta = pik_beta.array() + (double)1/n_k;


	MatrixXd post_gamma_beta = gamma_beta;
	gamma_beta = gamma_beta.array() + (double)1/n_k;

	VectorXd sigma2k (n_k);
	VectorXd vk (n_k);
	VectorXd Elogsigmak (n_k);
	VectorXd Elogvk (n_k);
	VectorXd sumElogvl (n_k);

	VectorXd pikexp1(n_k);
	VectorXd pikexp2(n_k);
	VectorXd index;
	VectorXd a_k = (VectorXd::Random(n_k)).array().abs();
	VectorXd b_k = (VectorXd::Random(n_k)).array().abs();
	VectorXd  kappa_k = pik_beta.colwise().sum();
	VectorXd lambda_k = (VectorXd::Random(n_k)).array().abs();

	XEbeta.setZero();//set to zero
	//for (size_t i=0; i<n_snp; i++)   {
	for (size_t i=0; i<T_snp; i++)   {
		//Eigen_matrix_get_row (UtX, i, x_col);
		//x_col   =  UtX.col(snp_no(i));
		//XEbeta += x_col*beta(snp_no(i));
		x_col   =  UtX.col(i);
		XEbeta += x_col*beta(i);
		}

	WEalpha.setZero(); //set to zero
	for (size_t j=0;j<n_j; j++)   {
		WEalpha += UtW.col(j)*m_alpha(j);
		}

	VectorXd y_res  = Uty - WEalpha - XEbeta;
	double sigma2e0 = ((y_res).dot(y_res))/(n_idv-1);
	cout << endl << endl;
	//double a_e  = (2*n_idv + n_snp)/2;
	double a_e  = (2*n_idv + T_snp)/2;
	double b_e  = (a_e-1)*sigma2e0;

	double ak   = 21;
	if (lambda < 0.01)     {lambda =    0.01;}
	else if (lambda > 100) {lambda =     100;}
	else                   {lambda =  lambda;}
	//double bk0  = lambda*(ak-1)/n_snp;
	double bk0  = lambda*(ak-1)/UtX.cols();

	VectorXd bk(n_k);
	bk(0) = bk0;
	for (size_t i=1;i<n_k; i++)  {
		bk(i) = bk(i-1)*2.7*sqrt(pow(i,i));
		}

	//for (size_t i=0;i<n_snp; i++)  {
	for (size_t i=0;i<T_snp; i++)  {
		for (size_t k=0;k<n_k; k++)  {
			sik2_beta(i,k) = bk(k)/(ak-1)*sigma2e0;
		}
		}

	double ae   = 0.1;
	double be   = 0.1;
	double a0   = 400;
	double b0   = 40;
	double a_lambda = a0 + n_k;
	double b_lambda = b0;
	double Elogsigmae, tx_ywx_res, xtxabk, A, post_Gn=0;
	double sigma2b = 0.001;

	double post_llike=0, post_sigma2e=0, post_sigma2b=0;

	VectorXd sigma2bX = (VectorXd::Random(w_step+s_step+1)).array().abs();
	sigma2bX(0)=sigma2b;

	VectorXd h = (VectorXd::Random(w_step+s_step+1)).array().abs();
	VectorXd H0(n_idv);
	VectorXd H(n_idv);

	VectorXd    llike = VectorXd::Zero(w_step+s_step);


	////initial values for a_k, b_k and sigma2k
	Ebeta2k = (mik_beta.cwiseProduct(mik_beta)) + sik2_beta;
	for (size_t k=0; k<n_k; k++)  {
		a_k(k) = (pik_beta.col(k).sum())/2 + ak;
		b_k(k) = (Ebeta2k.col(k).sum())*(a_e/b_e)/2 + bk(k);
		sigma2k(k) = b_k(k)/(a_k(k)-1);
		//cout<<sigma2k(k)<<endl;
		}

	A = (y_res).dot(y_res);
	//B1=(1/sigma2k.array()).rowwise().replicate(n_snp);
	B1=(1/sigma2k.array()).rowwise().replicate(T_snp);
	
	B1.col(0).fill(0);

	Ebeta2k = Ebeta2k.cwiseProduct(B1.transpose());
	double B = Ebeta2k.sum();
	a_e = n_idv + T_snp*0.01 + ae;
	b_e = (A + B + 2*be)/2;
	////random seed
	gsl_rng * rs = gsl_rng_alloc (gsl_rng_mt19937);

	if (n_k == 1)   {
		WritelmmBeta(beta);
		pheno_mean = m_alpha(0);
		}

	// //  begin MCMC sampling
	else {
	for (size_t S=0; S<(w_step + s_step); S++)   { 

		sigma2e    = 1/gsl_ran_gamma(rs,a_e,1/b_e);
		Elogsigmae = log(sqrt(sigma2e));

		// save Ebeta for the mixture normal component
		Ebeta = (beta_beta.cwiseProduct(gamma_beta)).rowwise().sum();
		
		if (S > (w_step-1))   {post_Ebeta += Ebeta;}

		//sample sigma2k and compute related quantities
		for (size_t k=0; k<n_k; k++)   {
			if (k==0) {
				sigma2k(k)    = 0; 
				Elogsigmak(k) = 0;
				}

			else {
				sigma2k(k)    = 1/gsl_ran_gamma(rs,a_k(k),1/b_k(k));
				Elogsigmak(k) = log(sqrt(sigma2k(k)));
				}

				if (k == (n_k-1)) {vk(k) = 0; Elogvk(k) = 0;}
				else {
					vk(k)     = gsl_ran_beta(rs,kappa_k(k),lambda_k(k));
					Elogvk(k) = log(vk(k));
				}
				sumElogvl(k)  = sum_Elogvl2(vk,k);
			}

		//////////////  sampling the mixture snp effects 
		XEbeta.setZero(); //set Gbeta to zero first
		//for (size_t i=0; i<n_snp; i++)   {
		//	x_col   =  UtX.col(snp_no(i));
		for (size_t i=0; i<T_snp; i++)   {
			x_col   =  UtX.col(i);
			XEbeta += x_col*Ebeta(i);
		}

		y_res.setZero();
		y_res = Uty - WEalpha;

		H0 = (sigma2b*D).array()+1;
		H0 = 1/H0.array();


		cout<<" S%xstep "<<S<<"  "<<S%50<<endl;
		//////////////////////////////////////
		//for (size_t i=0; i<n_snp; i++)   {
		for (size_t i=0; i<T_snp; i++)   {
			//x_col      = UtX.col(snp_no(i));
			x_col      = UtX.col(i);
			XEbeta    -= x_col*Ebeta(i);
			H = x_col.array()*H0.array();
			tx_ywx_res = H.dot(y_res - XEbeta);
			for (size_t k=0; k<n_k; k++)   {
				if (k == 0)   {
					mik_beta(i,k)  = 0; 
					sik2_beta(i,k) = 0; 
					pikexp1(k)     = 0;
					}
				else {
					xtxabk         = H.dot(x_col) + 1/sigma2k(k);
					mik_beta(i,k)  = tx_ywx_res/xtxabk;
					sik2_beta(i,k) = sigma2e/xtxabk;
					pikexp1(k)     = mik_beta(i,k)*mik_beta(i,k)/(2*sik2_beta(i,k)) + log(sqrt(sik2_beta(i,k))) - Elogsigmae - Elogsigmak(k);
				}
				beta_beta(i,k) = gsl_ran_gaussian(rs,sqrt(sik2_beta(i,k))) + mik_beta(i,k);
				pikexp2(k)     = Elogvk(k) + sumElogvl(k);
				}
			index = pikexp1 + pikexp2;
			index = (index.array() - index.maxCoeff()).exp();
			pik_beta.row(i) = index/index.sum();
		// multinomial sampling
		double mult_prob[n_k];
		unsigned int mult_no[n_k];
		for (size_t k=0; k<n_k; k++)   {
			mult_prob[k] = pik_beta(i,k);
			}
		gsl_ran_multinomial(rs, n_k, 1, mult_prob, mult_no);
		for (size_t k=0; k<n_k; k++)   {
			gamma_beta(i,k) = mult_no[k];
		}
		// multinomial sampling
		Ebeta(i) = ((beta_beta.row(i)).dot(gamma_beta.row(i)));
		XEbeta  += x_col*Ebeta(i);
		}
		//////////////////////////////////////




		if (S > (w_step-1))   {
			post_gamma_beta += gamma_beta;
			}
		//////////////  sampling the mixture snp effects 

		WEalpha.setZero();
		for (size_t j=0;j<n_j; j++)   {
			WEalpha += UtW.col(j)*Ealpha(j);
		}
		
		y_res.setZero();
		y_res = Uty - XEbeta;
		for (size_t j=0;j<n_j; j++)   {
			WEalpha    -= UtW.col(j)*Ealpha(j);
					  H = UtW.col(j).array()*H0.array();
			m_alpha(j)  = (H.dot(y_res - WEalpha))/(H.dot(UtW.col(j)));
			s2_alpha(j) = sigma2e/(H.dot(UtW.col(j)));
			Ealpha(j)   = m_alpha(j) + gsl_ran_gaussian(rs,sqrt(s2_alpha(j)));
			WEalpha    += UtW.col(j)*Ealpha(j);
			}

		if (S>(w_step-1))  {post_Ealpha += Ealpha;}

		a_lambda = a0 + n_k;
		b_lambda = b0 - sum_Elogvl2(vk,n_k-1);
		lambdax  = gsl_ran_gamma(rs,a_lambda,1/b_lambda);

		for (size_t k=0; k<n_k-1; k++)   {
			//kappa_k(k) = gamma_beta.col(k).sum()+(T_snp-n_snp)*(k==0)+1;
			kappa_k(k) = gamma_beta.col(k).sum()+1;
			lambda_k(k) = sum_lambda_k(gamma_beta,k,n_k) + lambdax;
			}

		Ebeta2k = beta_beta.cwiseProduct(beta_beta).cwiseProduct(gamma_beta);
		for (size_t k=0; k<n_k; k++)   {
			//a_k(k) = (gamma_beta.col(k).sum()+(T_snp-n_snp)*(k==0))/2 + ak;
			a_k(k) = (gamma_beta.col(k).sum())/2 + ak;
			b_k(k) = (Ebeta2k.col(k).sum())/(2*sigma2e) + bk(k);
			}

		y_res.setZero();
		y_res = Uty - XEbeta - WEalpha;

		bv0 = (y_res*sigma2b).array()*H0.array();
		if (S>(w_step-1))   {bv += bv0;}

		H = y_res.array()*H0.array();

		A  = H.dot(y_res);
		B1 = (1/sigma2k.tail(n_k-1).array()).rowwise().replicate(T_snp);
		//B1 = (1/sigma2k.tail(n_k-1).array()).rowwise().replicate(n_snp);
		Ebeta2k   = Ebeta2k.rightCols(n_k-1).cwiseProduct(B1.transpose());
		double B  = Ebeta2k.sum();
		double Gn = gamma_beta.rightCols(n_k-1).sum();
		
		a_e = n_idv/2 + Gn/2 + ae;
		b_e = (A + B + 2*be)/2;

		//double h_new  = gsl_rng_uniform(rs);
		double h_new  = gsl_ran_beta(rs,2,8);
		double ratio1 = log_h(D,y_res,h_new,sigma2e,ae,be)-log(gsl_ran_beta_pdf(h_new,2,8));
		double ratio2 = log_h(D,y_res,h(S), sigma2e,ae,be)-log(gsl_ran_beta_pdf(h(S),2,8));
		double ratio  = exp(ratio1-ratio2);
		if (ratio>1) {ratio = 1;}
		else         {ratio = ratio;}
		double u = gsl_rng_uniform(rs);
		if (u<ratio) {h(S+1) = h_new;}
		else         {h(S+1) = h(S);}

		sigma2b = h(S+1)/(1-h(S+1));

		if (S > (w_step-1)) {post_Gn += Gn;}
		////////////////////////////////////
		double llike0=logLike(D,y_res,sigma2b,sigma2e); 
		     llike(S)=llike0; 

		if (S > (w_step-1))   {	
			post_llike   += llike0;
			post_sigma2e += sigma2e;
			post_sigma2b += sigma2b;
			}

		size_t xstep = (s_step + w_step)*0.2+1;
		if (S%xstep==0||S==(s_step + w_step-1)) {ProgressBar ("screen MCMC sampling ", S, s_step+w_step-1);}
		}
	}
	
	mixture_no = snp_no.size();
	cout<<endl<<endl<<"MCMC sampling is finished and now compute the polygenic effects ..."<<endl;

	VectorXd eigen_alpha = VectorXd::Zero (T_snp);
	VectorXd Beta        = VectorXd::Zero (T_snp);
	eigen_alpha          = UtX.transpose()*(bv/s_step)/T_snp;

	//int mcount=0;
	//for (size_t i=0; i<T_snp; i++)  {
	//	if (mcount<snp_no.size())  {
	//		if (i==(size_t)snp_no(mcount))  {
	//			Beta(i) = post_Ebeta(mcount)/s_step;
	//			mcount++;
	//			}}
	//	else {Beta(i) = Beta(i);}
	//}

	WriteCoeff (eigen_alpha, post_Ebeta/s_step);
	pheno_mean = post_Ealpha(0)/s_step;

/////////
	XEbeta.setZero(); //set Gbeta to zero first
		for (size_t i=0; i<T_snp; i++)   {
		//for (size_t i=0; i<n_snp; i++)   {
			x_col   =  UtX.col(i);
			XEbeta += x_col*post_Ebeta(i)/s_step;
		}
	WEalpha.setZero();
		for (size_t j=0;j<n_j; j++)   {
			WEalpha += UtW.col(j)*post_Ealpha(j)/s_step;
		}

	y_res.setZero();
	y_res = Uty - XEbeta - WEalpha;

	double llike_hat = logLike(D,y_res,post_sigma2b/s_step,post_sigma2e/s_step); 
	VectorXd llike2  = llike.tail(s_step).array();
	          llike2 = llike2.array() - llike2.sum()/s_step;
	pD1   =  2*(llike_hat-post_llike/s_step);
	pD2   =  2*(llike2.dot(llike2))/(s_step-1);
	if (pD1 < 0)   {pD1 = 1;}
	DIC1  = -2* llike_hat +  2*pD1;
	DIC2  = -2* llike_hat +  2*pD2;
	BIC1  = -2* llike_hat + log(n_idv)*pD1;
	BIC2  = -2* llike_hat + log(n_idv)*pD2;

	cout<<" K     = " << n_k<<endl;
	cout<<" pD1   = " << pD1 <<endl;
	cout<<" pD2   = " << pD2<<endl;
	cout<<" DIC1  = " << DIC1<<endl;
	cout<<" DIC2  = " << DIC2<<endl;
	cout<<" BIC1  = " << BIC1<<endl;
	cout<<" BIC2  = " << BIC2<<endl;

	cout<<endl<< "Computaion Time for LDR.Gibbs = ";
	cout<<(clock()-time_begin)/(double(CLOCKS_PER_SEC)*60.0)<<" min"<<endl<<endl;

	////////////
	xtx.resize(0);xty.resize(0);wtw.resize(0);wty.resize(0);Utu.resize(0);
	Ue.resize(0);Ub.resize(0);bv0.resize(0);bv.resize(0);V.resize(0);
	M.resize(0);Ealpha.resize(0);post_Ealpha.resize(0);m_alpha.resize(0);
	s2_alpha.resize(0);WEalpha.resize(0);XEbeta.resize(0);Ebeta.resize(0);
	post_Ebeta.resize(0);Ebeta2k.resize(0,0);B1.resize(0,0);mik_beta.resize(0,0);
	sik2_beta.resize(0,0);pik_beta.resize(0,0);beta_beta.resize(0,0);
	gamma_beta.resize(0,0);sigma2k.resize(0);vk.resize(0);Elogsigmak.resize(0);
	Elogvk.resize(0);sumElogvl.resize(0);pikexp1.resize(0);pikexp2.resize(0);
	index.resize(0);a_k.resize(0);b_k.resize(0);kappa_k.resize(0);
	lambda_k.resize(0);y_res.resize(0);bk.resize(0);Beta.resize(0);
}


void LDR::Gibbs_without_u_screen (const MatrixXd UtX, const VectorXd Uty, const MatrixXd UtW,const VectorXd D, const VectorXd Wbeta, const VectorXd se_Wbeta, const VectorXd beta, const VectorXi snp_no, double lambda) 
{
	
	setprecision(5); 
	clock_t time_begin  = clock();
	size_t n_snp = UtX.cols(); //n * p
	size_t s_snp = snp_no.size();
	size_t n_idv = Uty.size();
	size_t n_j   = UtW.cols();
	VectorXd x_col(n_idv);

	VectorXd xtx(n_snp);
	VectorXd xty(n_snp);
	VectorXd wtw(n_j);
	VectorXd wty(n_j);
	VectorXd Utu = VectorXd::Zero(n_idv);
	VectorXd Ue(n_idv);
	VectorXd Ub(n_idv);

	VectorXd bv0 = VectorXd::Zero(n_idv); //for computation of alpha
	VectorXd bv  = VectorXd::Zero(n_idv); //for computation of alpha
	VectorXd V = VectorXd::Zero (n_idv);//compute Utu 
	VectorXd M = VectorXd::Zero (n_idv);//compute Utu
	VectorXd snp_label(n_snp);

	VectorXd Ealpha      = Wbeta;//intercept
	VectorXd post_Ealpha = VectorXd::Zero(n_j);//save intercept
	VectorXd m_alpha     = Wbeta;
	VectorXd s2_alpha    = se_Wbeta.array()*se_Wbeta.array();

	VectorXd WEalpha(n_idv); 
	double   sigma2e=0, lambdax;
	VectorXd XEbeta(n_idv); //G*beta
	VectorXd Ebeta(n_snp);
	VectorXd post_Ebeta = VectorXd::Zero(n_snp);//save beta
	MatrixXd Ebeta2k;
	MatrixXd B1;
	MatrixXd mik_beta = MatrixXd::Zero(n_snp, n_k);

	VectorXd    llike = VectorXd::Zero(w_step+s_step);

	for (size_t i=0;i<n_snp; i++)  {
		//Eigen_matrix_get_row (UtX, i, x_col);
		x_col  = UtX.col(i);
		xtx(i) = x_col.dot(x_col);
		xty(i) = x_col.dot(Uty);
	}

	for (size_t j=0;j<n_j; j++)  {
		wtw(j) = UtW.col(j).dot(UtW.col(j));
		wty(j) = UtW.col(j).dot(Uty);
	}

	//// initial values for SNP and each normal
	//// component has the same initial values
	for (size_t i=0;i<n_snp; i++)  {
		for (size_t k=1;k<n_k; k++)  {
			mik_beta(i,k) = beta(i);
		}
	}

	mik_beta.col(0).fill(0);
	MatrixXd sik2_beta  =(MatrixXd::Random(n_snp, n_k)).array().abs();
	MatrixXd pik_beta   = MatrixXd::Zero  (n_snp, n_k);
	pik_beta = pik_beta.array() + (double)1/n_k;

	MatrixXd beta_beta  = MatrixXd::Zero(n_snp, n_k);
	MatrixXd gamma_beta = MatrixXd::Zero(n_snp, n_k);
	MatrixXd post_gamma_beta = gamma_beta;
	gamma_beta = gamma_beta.array() + (double)1/n_k;

	VectorXd sigma2k (n_k);
	VectorXd vk (n_k);
	VectorXd Elogsigmak (n_k);
	VectorXd Elogvk (n_k);
	VectorXd sumElogvl (n_k);

	VectorXd pikexp1(n_k);
	VectorXd pikexp2(n_k);
	VectorXd index;
	VectorXd a_k = (VectorXd::Random(n_k)).array().abs();
	VectorXd b_k = (VectorXd::Random(n_k)).array().abs();
	VectorXd  kappa_k = pik_beta.colwise().sum();
	VectorXd lambda_k = (VectorXd::Random(n_k)).array().abs();

	XEbeta.setZero();//set to zero
	//for (size_t i=0; i<n_snp; i++)   {
		//x_col   =  UtX.col(i);
		//XEbeta += x_col*beta(i);
		//}
	XEbeta = UtX*Ebeta;

	WEalpha.setZero(); //set to zero
	for (size_t j=0;j<n_j; j++)   {
		WEalpha += UtW.col(j)*m_alpha(j);
		}

	VectorXd y_res  = Uty - WEalpha - XEbeta;
	double sigma2e0 = ((y_res).dot(y_res))/(n_idv-1);
	double a_e  = (2*n_idv + n_snp)/2;
	double b_e  = (a_e-1)*sigma2e0;

	double ak   = 21;
	if (lambda < 0.01)     {lambda =    0.01;}
	else if (lambda > 100) {lambda =     100;}
	else                   {lambda =  lambda;}
	double bk0  = lambda*(ak-1)/n_snp;

	VectorXd bk(n_k);
	bk(0) = bk0;
	for (size_t i=1;i<n_k; i++)  {
		bk(i) = bk(i-1)*1.7*sqrt(pow(i,i));
		}

	for (size_t i=0;i<n_snp; i++)  {
		for (size_t k=0;k<n_k; k++)  {
			sik2_beta(i,k) = bk(k)/(ak-1)*sigma2e0;
		}
	}

	double ae   = 0.1;
	double be   = 0.1;
	//double a0   = 100;
	//double b0   = 10;
	double a0   = 1;
	double b0   = 0.1;
	double a_lambda = a0 + n_k;
	double b_lambda = b0;
	double Elogsigmae, tx_ywx_res, xtxabk, A, post_Gn=0;
	double sigma2b = 0.1;

	double post_llike=0, post_sigma2e=0, post_sigma2b=0;

	VectorXd sigma2bX = (VectorXd::Random(w_step+s_step+1)).array().abs();
	sigma2bX(0)=sigma2b;

	VectorXd h = (VectorXd::Random(w_step+s_step+1)).array().abs();
	VectorXd H0(n_idv);
	VectorXd H(n_idv);

	//for (size_t k=0; k<(w_step+s_step+1); k++)  {
		//cout<<h(k)<<endl;
	//}

	////initial values for a_k, b_k and sigma2k
	Ebeta2k = (mik_beta.cwiseProduct(mik_beta)) + sik2_beta;
	for (size_t k=0; k<n_k; k++)  {
		a_k(k) = (pik_beta.col(k).sum())/2 + ak;
		b_k(k) = (Ebeta2k.col(k).sum())*(a_e/b_e)/2 + bk(k);
		sigma2k(k) = b_k(k)/(a_k(k)-1);
		}

	A = (y_res).dot(y_res);
	B1=(1/sigma2k.array()).rowwise().replicate(n_snp);
	B1.col(0).fill(0);

	Ebeta2k = Ebeta2k.cwiseProduct(B1.transpose());
	double B = Ebeta2k.sum();
	a_e = n_idv + n_snp*0.1 + ae;
	b_e = (A + B + 2*be)/2;

	////random seed
	gsl_rng * rs = gsl_rng_alloc (gsl_rng_mt19937);

	if (n_k == 1)   {
		WritelmmBeta(beta);
		pheno_mean = m_alpha(0);
		}

	// //  begin MCMC sampling
	else {
	for (size_t S=0; S<(w_step + s_step); S++)   {

		sigma2e    = 1/gsl_ran_gamma(rs,a_e,1/b_e);
		Elogsigmae = log(sqrt(sigma2e));

		// save Ebeta for the mixture normal component
		Ebeta = (beta_beta.cwiseProduct(gamma_beta)).rowwise().sum();

		if (S > (w_step-1))   {post_Ebeta += Ebeta;}

		//sample sigma2k and compute related quantities
		for (size_t k=0; k<n_k; k++)   {
			if (k==0) {
				sigma2k(k)    = 0; 
				Elogsigmak(k) = 0;
				}

			else {
				sigma2k(k)    = 1/gsl_ran_gamma(rs,a_k(k),1/b_k(k));
				Elogsigmak(k) = log(sqrt(sigma2k(k)));
				}

				if (k == (n_k-1)) {vk(k) = 0; Elogvk(k) = 0;}
				else {
					vk(k)     = gsl_ran_beta(rs,kappa_k(k),lambda_k(k));
					Elogvk(k) = log(vk(k));
				}
				sumElogvl(k)  = sum_Elogvl2(vk,k);
			}

		//////////////  sampling the mixture snp effects 
		XEbeta.setZero(); //set Gbeta to zero first
			//for (size_t i=0; i<n_snp; i++)   {
				//x_col   =  UtX.col(i);
				//XEbeta += x_col*Ebeta(i);
		//}
		XEbeta = UtX*Ebeta;

		y_res.setZero();
		y_res = Uty - WEalpha;

		H0 = (sigma2b*D).array()+1;
		H0 = 1/H0.array();

		//////////////////////////////////////
		// full sampling, i.e., sampling effects for all the snps 
		//cout<<tp<<endl;
		size_t tpx=tp;
		if (S%tpx==0)  {
		for (size_t i=0; i<n_snp; i++)   {
			x_col      = UtX.col(i);
			XEbeta    -= x_col*Ebeta(i);

			H = x_col.array()*H0.array();
			tx_ywx_res = H.dot(y_res - XEbeta);

			for (size_t k=0; k<n_k; k++)   {

				if (k == 0) {
					mik_beta(i,k)  = 0; 
					sik2_beta(i,k) = 0; 
					pikexp1(k)     = 0;
					}

				else {
					xtxabk         = H.dot(x_col) + 1/sigma2k(k);
					mik_beta(i,k)  = tx_ywx_res/xtxabk;
					sik2_beta(i,k) = sigma2e/xtxabk;
					pikexp1(k)     = mik_beta(i,k)*mik_beta(i,k)/(2*sik2_beta(i,k)) + log(sqrt(sik2_beta(i,k))) - Elogsigmae - Elogsigmak(k);
				}

				beta_beta(i,k) = gsl_ran_gaussian(rs,sqrt(sik2_beta(i,k))) + mik_beta(i,k);
				pikexp2(k)     = Elogvk(k) + sumElogvl(k);
				}

			index = pikexp1 + pikexp2;
			index = (index.array() - index.maxCoeff()).exp();
			pik_beta.row(i) = index/index.sum();

		// multinomial sampling
		double mult_prob[n_k];
		unsigned int mult_no[n_k];

		for (size_t k=0; k<n_k; k++)   {
			mult_prob[k] = pik_beta(i,k);
			}

		gsl_ran_multinomial(rs, n_k, 1, mult_prob, mult_no);
		for (size_t k=0; k<n_k; k++)   {
			gamma_beta(i,k) = mult_no[k];
		}

		Ebeta(i) = ((beta_beta.row(i)).dot(gamma_beta.row(i)));
		XEbeta  += x_col*Ebeta(i);
		}
		}

		else {
			//partial sampling
			//i.e., sampling effects for pre-selected snps with large effects
			for (size_t t=0; t<s_snp; t++)   {
			  size_t i = snp_no(t);
			x_col      = UtX.col(i);
			XEbeta    -= x_col*Ebeta(i);

			H = x_col.array()*H0.array();
			tx_ywx_res = H.dot(y_res - XEbeta);

			for (size_t k=0; k<n_k; k++)   {

				if (k == 0) {
					mik_beta(i,k)  = 0; 
					sik2_beta(i,k) = 0; 
					pikexp1(k)     = 0;
					}

				else {
					xtxabk         = H.dot(x_col) + 1/sigma2k(k);
					mik_beta(i,k)  = tx_ywx_res/xtxabk;
					sik2_beta(i,k) = sigma2e/xtxabk;
					pikexp1(k)     = mik_beta(i,k)*mik_beta(i,k)/(2*sik2_beta(i,k)) + log(sqrt(sik2_beta(i,k))) - Elogsigmae - Elogsigmak(k);
				}

				beta_beta(i,k) = gsl_ran_gaussian(rs,sqrt(sik2_beta(i,k))) + mik_beta(i,k);
				pikexp2(k)     = Elogvk(k) + sumElogvl(k);
				}

			index = pikexp1 + pikexp2;
			index = (index.array() - index.maxCoeff()).exp();
			pik_beta.row(i) = index/index.sum();

		// multinomial sampling
		double mult_prob[n_k];
		unsigned int mult_no[n_k];

		for (size_t k=0; k<n_k; k++)   {
			mult_prob[k] = pik_beta(i,k);
			}

		gsl_ran_multinomial(rs, n_k, 1, mult_prob, mult_no);
		for (size_t k=0; k<n_k; k++)   {
			gamma_beta(i,k) = mult_no[k];
		}

		Ebeta(i) = ((beta_beta.row(i)).dot(gamma_beta.row(i)));
		XEbeta  += x_col*Ebeta(i);
		}

		}
		//////////////////////////////////////


		if (S > (w_step-1))   {post_gamma_beta += gamma_beta;}
		//////////////  sampling the mixture snp effects 

		WEalpha.setZero();
		for (size_t j=0;j<n_j; j++)   {
			WEalpha += UtW.col(j)*Ealpha(j);
		}

		y_res.setZero();
		y_res = Uty - XEbeta;
		for (size_t j=0;j<n_j; j++)   {
			WEalpha    -= UtW.col(j)*Ealpha(j);
					  H = UtW.col(j).array()*H0.array();
			m_alpha(j)  = (H.dot(y_res - WEalpha))/(H.dot(UtW.col(j)));
			s2_alpha(j) = sigma2e/(H.dot(UtW.col(j)));
			Ealpha(j)   = m_alpha(j) + gsl_ran_gaussian(rs,sqrt(s2_alpha(j)));
			WEalpha    += UtW.col(j)*Ealpha(j);
			}

		if (S>(w_step-1))   {post_Ealpha += Ealpha;}

		a_lambda = a0 + n_k;
		b_lambda = b0 - sum_Elogvl2(vk,n_k-1);
		lambdax  = gsl_ran_gamma(rs,a_lambda,1/b_lambda);

		for (size_t k=0; k<n_k-1; k++)   {
			kappa_k(k)  = gamma_beta.col(k).sum() + 1;
			lambda_k(k) = sum_lambda_k(gamma_beta,k,n_k) + lambdax;
			}

		Ebeta2k = beta_beta.cwiseProduct(beta_beta).cwiseProduct(gamma_beta);
		for (size_t k=0; k<n_k; k++)   {
			a_k(k) = (gamma_beta.col(k).sum())/2 + ak;
			b_k(k) = (Ebeta2k.col(k).sum())/(2*sigma2e) + bk(k);
			}

		/*y_res.setZero();
		y_res = Uty - XEbeta - WEalpha;
		for (size_t i=0; i<n_idv; i++) {
			V(i)  = sigma2b*D(i)/(sigma2b*D(i) + 1);
			M(i)  = y_res(i)*V(i);
			Utu(i)= M(i) + gsl_ran_gaussian(rs,sqrt(V(i)*sigma2e));
			if (D(i)  == 0) {
				Ue(i) = 0; 
				Ub(i) = 0;
				}
			else {
				Ue(i) = Utu(i)/(sigma2b*D(i)); //for sigma2e
				Ub(i) = Utu(i)/(sigma2e*D(i)); //for sigma2b
			}
			bv0(i) = y_res(i)*sigma2b/(sigma2b*D(i) + 1); 
			}

		if (S>(w_step-1))   {bv += bv0;}
		*/

		y_res.setZero();
		y_res = Uty - XEbeta - WEalpha;

		bv0 = (y_res*sigma2b).array()*H0.array();
		if (S>(w_step-1))   {bv += bv0;}

		H = y_res.array()*H0.array();

		A  = H.dot(y_res);
		B1 = (1/sigma2k.tail(n_k-1).array()).rowwise().replicate(n_snp);
		Ebeta2k   = Ebeta2k.rightCols(n_k-1).cwiseProduct(B1.transpose());
		double B  = Ebeta2k.sum();
		double Gn = gamma_beta.rightCols(n_k-1).sum();
		
		a_e = n_idv/2 + Gn/2 + ae;
		b_e = (A + B + 2*be)/2;

		/*double sigma2bX_new = sigma2bX(S) + gsl_ran_gaussian(rs,sqrt(1));
		double ratio1 = log_sigma2b(D,y_res,sigma2bX_new,sigma2e,ae,be);
		double ratio2 = log_sigma2b(D,y_res,sigma2bX(S), sigma2e,ae,be);
		double ratio  = exp(ratio1-ratio2);
		if (ratio>1) {ratio = 1;}
		else         {ratio = ratio;}
		if (sigma2bX_new<0) {ratio = 0;}
		*/

		//double h_new  = gsl_rng_uniform(rs);
		double h_new  = gsl_ran_beta(rs,2,8);
		double ratio1 = log_h(D,y_res,h_new,sigma2e,ae,be)-log(gsl_ran_beta_pdf(h_new,2,8));
		double ratio2 = log_h(D,y_res,h(S), sigma2e,ae,be)-log(gsl_ran_beta_pdf(h(S),2,8));
		double ratio  = exp(ratio1-ratio2);
		if (ratio>1) {ratio = 1;}
		else         {ratio = ratio;}
		double u = gsl_rng_uniform(rs);
		if (u<ratio) {h(S+1) = h_new;}
		else         {h(S+1) = h(S);}

		sigma2b = h(S+1)/(1-h(S+1));

		//cout<<S<<" "<<setprecision(6)<<h_new<<" "<<sigma2b<<" "<<sigma2e<<" "<<sigma2b*sigma2e<<endl;
		//cout<<setprecision(6)<<sigma2bX_new<<" "<<sigma2b<<" "<<sigma2e<<" "<<sigma2b*sigma2e<<endl;
		if (S > (w_step-1))   {post_Gn += Gn;}
		////////////////////////////////////

		double llike0=logLike(D,y_res,sigma2b,sigma2e); 
		     llike(S)=llike0; 

		if (S > (w_step-1))   {	
			post_llike   += llike0;
			post_sigma2e += sigma2e;
			post_sigma2b += sigma2b;
			}

		size_t xstep = (s_step + w_step)*0.2+1;
		if (S%xstep==0||S==(s_step + w_step-1)) {
			ProgressBar ("MCMC sampling ", S, s_step+w_step-1);
			}
		}
	}

	mixture_no = post_Gn/s_step;
	cout<<endl<<"MCMC sampling is finished"<<endl;
	VectorXd eigen_alpha = VectorXd::Zero (n_snp);
	eigen_alpha = UtX.transpose()*(bv/s_step)/n_snp;
	WriteCoeff (eigen_alpha, post_Ebeta/s_step);
	pheno_mean = post_Ealpha(0)/s_step;

/////////
	XEbeta.setZero(); //set Gbeta to zero first
		for (size_t i=0; i<n_snp; i++)   {
			x_col   =  UtX.col(i);
			XEbeta += x_col*post_Ebeta(i)/s_step;
		}
	WEalpha.setZero();
		for (size_t j=0;j<n_j; j++)   {
			WEalpha += UtW.col(j)*post_Ealpha(j)/s_step;
		}

	y_res.setZero();
	y_res = Uty - XEbeta - WEalpha;

	double llike_hat = logLike(D,y_res,post_sigma2b/s_step,post_sigma2e/s_step); 
	VectorXd llike2  = llike.tail(s_step).array();
	          llike2 = llike2.array() - llike2.sum()/s_step;
	pD1   =  2*(llike_hat-post_llike/s_step);
	pD2   =  2*(llike2.dot(llike2))/(s_step-1);
	if (pD1 < 0)   {pD1 = 1;}
	DIC1  = -2* llike_hat +  2*pD1;
	DIC2  = -2* llike_hat +  2*pD2;
	BIC1  = -2* llike_hat + log(n_idv)*pD1;
	BIC2  = -2* llike_hat + log(n_idv)*pD2;

	//cout<<" K     = " << n_k<<endl;
	//cout<<" pD1   = " << pD1 <<endl;
	//cout<<" pD2   = " << pD2<<endl;
	//cout<<" DIC1  = " << DIC1<<endl;
	//cout<<" DIC2  = " << DIC2<<endl;
	//cout<<" BIC1  = " << BIC1<<endl;
	//cout<<" BIC2  = " << BIC2<<endl;

	//ofstream outfile1;
	//outfile1.open ("llike.txt", ios::out | ios::binary);
	//for (size_t k = 0; k<s_step; k++)  {
	//	outfile1 << setprecision(7) <<llike2(k)<<endl;
	//	}
	//outfile1 << endl;

	//for (size_t k = 0; k<(w_step+s_step); k++)  {
	//	outfile1 << setprecision(7) <<llike(k)<<endl;
	//	}
	//outfile1 << endl;

	WriteLike(llike,w_step+s_step);

	cout<<"Computaion Time for DPR.Gibbs = ";
	cout<<(clock()-time_begin)/(double(CLOCKS_PER_SEC)*60.0)<<" min"<<endl<<endl;
/*
	//output mik and pik
	ofstream outfile1;
	ofstream outfile2;
	ofstream outfile3;
	ofstream outfile4;
	outfile1.open ("mik.txt", ios::out | ios::binary);
	outfile2.open ("pik.txt", ios::out | ios::binary);
	outfile3.open ("sik.txt", ios::out | ios::binary);
	//outfile4.open ("sigmak.txt", ios::out | ios::binary);
	
	for (size_t i = 0; i<n_snp; i++)  {
		for (size_t k = 0; k<n_k; k++)  {
			outfile1 << setprecision(4) <<  mik_beta(i,k)<<"\t";
			outfile2 << setprecision(4) <<  pik_beta(i,k)<<"\t";
			outfile3 << setprecision(4) << sik2_beta(i,k)<<"\t";
			//outfile4 << setprecision(4) << b_k(k)/a_k(k) <<"\t";
			}
		outfile1 << endl;
		outfile2 << endl;
		outfile3 << endl;
		//outfile4 << endl;
	}
	outfile1.close();
	outfile2.close();
	outfile3.close();
	//outfile4.close();
	//output mik and pik
*/

	////////////
	xtx.resize(0);xty.resize(0);wtw.resize(0);wty.resize(0);
	Utu.resize(0);Ue.resize(0);Ub.resize(0);bv0.resize(0);
	bv.resize(0);V.resize(0);M.resize(0);Ealpha.resize(0);
	post_Ealpha.resize(0);m_alpha.resize(0);s2_alpha.resize(0);
	WEalpha.resize(0);XEbeta.resize(0);Ebeta.resize(0);
	post_Ebeta.resize(0);Ebeta2k.resize(0,0);
	B1.resize(0,0);mik_beta.resize(0,0);sik2_beta.resize(0,0);
	pik_beta.resize(0,0);beta_beta.resize(0,0);gamma_beta.resize(0,0);
	sigma2k.resize(0);vk.resize(0);Elogsigmak.resize(0);
	Elogvk.resize(0);sumElogvl.resize(0);pikexp1.resize(0);
	pikexp2.resize(0);index.resize(0);a_k.resize(0);
	b_k.resize(0);kappa_k.resize(0);lambda_k.resize(0);
	y_res.resize(0);bk.resize(0);
	eigen_alpha.resize(0);
}



void LDR::Gibbs_without_u_screen_dic0 (const MatrixXd UtX, const VectorXd Uty, const MatrixXd UtW,const VectorXd D, const VectorXd Wbeta, const VectorXd se_Wbeta, const VectorXd beta, const VectorXi snp_no, double lambda,size_t n_k) 
{
	
	setprecision(5); 
	size_t w_step0=w_step*sp;
	size_t s_step0=s_step*sp;
	//clock_t time_begin  = clock();
	size_t n_snp = UtX.cols(); //n * p
	size_t s_snp = snp_no.size();
	size_t n_idv = Uty.size();
	size_t n_j   = UtW.cols();
	VectorXd x_col(n_idv);

	VectorXd xtx(n_snp);
	VectorXd xty(n_snp);
	VectorXd wtw(n_j);
	VectorXd wty(n_j);
	VectorXd Utu = VectorXd::Zero(n_idv);
	VectorXd Ue(n_idv);
	VectorXd Ub(n_idv);

	VectorXd bv0 = VectorXd::Zero(n_idv); //for computation of alpha
	VectorXd bv  = VectorXd::Zero(n_idv); //for computation of alpha
	VectorXd V = VectorXd::Zero (n_idv);//compute Utu 
	VectorXd M = VectorXd::Zero (n_idv);//compute Utu
	VectorXd snp_label(n_snp);

	VectorXd Ealpha      = Wbeta;//intercept
	VectorXd post_Ealpha = VectorXd::Zero(n_j);//save intercept
	VectorXd m_alpha     = Wbeta;
	VectorXd s2_alpha    = se_Wbeta.array()*se_Wbeta.array();

	VectorXd WEalpha(n_idv); 
	double   sigma2e=0, lambdax;
	VectorXd XEbeta(n_idv); //G*beta
	VectorXd Ebeta(n_snp);
	VectorXd post_Ebeta = VectorXd::Zero(n_snp);//save beta
	MatrixXd Ebeta2k;
	MatrixXd B1;
	MatrixXd mik_beta = MatrixXd::Zero(n_snp, n_k);
	VectorXd    llike = VectorXd::Zero(w_step0+s_step0);

	for (size_t i=0;i<n_snp; i++)  {
		//Eigen_matrix_get_row (UtX, i, x_col);
		x_col  = UtX.col(i);
		xtx(i) = x_col.dot(x_col);
		xty(i) = x_col.dot(Uty);
	}

	for (size_t j=0;j<n_j; j++)  {
		wtw(j) = UtW.col(j).dot(UtW.col(j));
		wty(j) = UtW.col(j).dot(Uty);
	}

	//// initial values for SNP and each normal
	//// component has the same initial values
	for (size_t i=0;i<n_snp; i++)  {
		for (size_t k=1;k<n_k; k++)  {
			mik_beta(i,k) = beta(i);
		}
	}

	mik_beta.col(0).fill(0);
	MatrixXd sik2_beta  =(MatrixXd::Random(n_snp, n_k)).array().abs();
	MatrixXd pik_beta   = MatrixXd::Zero  (n_snp, n_k);
	pik_beta = pik_beta.array() + (double)1/n_k;

	MatrixXd beta_beta  = MatrixXd::Zero(n_snp, n_k);
	MatrixXd gamma_beta = MatrixXd::Zero(n_snp, n_k);
	MatrixXd post_gamma_beta = gamma_beta;
	gamma_beta = gamma_beta.array() + (double)1/n_k;

	VectorXd sigma2k (n_k);
	VectorXd vk (n_k);
	VectorXd Elogsigmak (n_k);
	VectorXd Elogvk (n_k);
	VectorXd sumElogvl (n_k);

	VectorXd pikexp1(n_k);
	VectorXd pikexp2(n_k);
	VectorXd index;
	VectorXd a_k = (VectorXd::Random(n_k)).array().abs();
	VectorXd b_k = (VectorXd::Random(n_k)).array().abs();
	VectorXd  kappa_k = pik_beta.colwise().sum();
	VectorXd lambda_k = (VectorXd::Random(n_k)).array().abs();

	XEbeta.setZero();//set to zero
	//for (size_t i=0; i<n_snp; i++)   {
		//x_col   =  UtX.col(i);
		//XEbeta += x_col*beta(i);
		//}
	XEbeta = UtX*Ebeta;

	WEalpha.setZero(); //set to zero
	for (size_t j=0;j<n_j; j++)   {
		WEalpha += UtW.col(j)*m_alpha(j);
		}

	VectorXd y_res  = Uty - WEalpha - XEbeta;
	double sigma2e0 = ((y_res).dot(y_res))/(n_idv-1);
	//cout << endl << endl;
	double a_e  = (2*n_idv + n_snp)/2;
	double b_e  = (a_e-1)*sigma2e0;

	double ak   = 21;
	if (lambda < 0.01)     {lambda =    0.01;}
	else if (lambda > 100) {lambda =     100;}
	else                   {lambda =  lambda;}
	double bk0  = lambda*(ak-1)/n_snp;

	VectorXd bk(n_k);
	bk(0) = bk0;
	for (size_t i=1;i<n_k; i++)  {
		bk(i) = bk(i-1)*1.7*sqrt(pow(i,i));
		}

	for (size_t i=0;i<n_snp; i++)  {
		for (size_t k=0;k<n_k; k++)  {
			sik2_beta(i,k) = bk(k)/(ak-1)*sigma2e0;
		}
	}

	double ae   = 0.1;
	double be   = 0.1;
	//double a0   = 40;
	//double b0   = 4;
	double a0   = 1;
	double b0   = 0.1;
	double a_lambda = a0 + n_k;
	double b_lambda = b0;
	double Elogsigmae, tx_ywx_res, xtxabk, A, post_Gn=0;
	double sigma2b = 0.1;

	double post_llike=0, post_sigma2e=0, post_sigma2b=0;

	VectorXd sigma2bX = (VectorXd::Random(w_step0+s_step0+1)).array().abs();
	sigma2bX(0)=sigma2b;

	VectorXd h = (VectorXd::Random(w_step0+s_step0+1)).array().abs();
	VectorXd H0(n_idv);
	VectorXd H(n_idv);

	//for (size_t k=0; k<(w_step0+s_step0+1); k++)  {
		//cout<<h(k)<<endl;
	//}

	////initial values for a_k, b_k and sigma2k
	Ebeta2k = (mik_beta.cwiseProduct(mik_beta)) + sik2_beta;
	for (size_t k=0; k<n_k; k++)  {
		a_k(k) = (pik_beta.col(k).sum())/2 + ak;
		b_k(k) = (Ebeta2k.col(k).sum())*(a_e/b_e)/2 + bk(k);
		sigma2k(k) = b_k(k)/(a_k(k)-1);
		}

	A = (y_res).dot(y_res);
	B1=(1/sigma2k.array()).rowwise().replicate(n_snp);
	B1.col(0).fill(0);

	Ebeta2k = Ebeta2k.cwiseProduct(B1.transpose());
	double B = Ebeta2k.sum();
	a_e = n_idv + n_snp*0.1 + ae;
	b_e = (A + B + 2*be)/2;

	////random seed
	gsl_rng * rs = gsl_rng_alloc (gsl_rng_mt19937);

	if (n_k == 1)   {
		WritelmmBeta(beta);
		pheno_mean = m_alpha(0);
		}

	// //  begin MCMC sampling
	else {
	for (size_t S=0; S<(w_step0 + s_step0); S++)   {

		sigma2e    = 1/gsl_ran_gamma(rs,a_e,1/b_e);
		Elogsigmae = log(sqrt(sigma2e));

		// save Ebeta for the mixture normal component
		Ebeta = (beta_beta.cwiseProduct(gamma_beta)).rowwise().sum();

		if (S > (w_step0-1))   {post_Ebeta += Ebeta;}

		//sample sigma2k and compute related quantities
		for (size_t k=0; k<n_k; k++)   {
			if (k==0) {
				sigma2k(k)    = 0; 
				Elogsigmak(k) = 0;
				}

			else {
				sigma2k(k)    = 1/gsl_ran_gamma(rs,a_k(k),1/b_k(k));
				Elogsigmak(k) = log(sqrt(sigma2k(k)));
				}

				if (k == (n_k-1)) {vk(k) = 0; Elogvk(k) = 0;}
				else {
					vk(k)     = gsl_ran_beta(rs,kappa_k(k),lambda_k(k));
					Elogvk(k) = log(vk(k));
				}
				sumElogvl(k)  = sum_Elogvl2(vk,k);
			}

		//////////////  sampling the mixture snp effects 
		XEbeta.setZero(); //set Gbeta to zero first
			//for (size_t i=0; i<n_snp; i++)   {
				//x_col   =  UtX.col(i);
				//XEbeta += x_col*Ebeta(i);
		//}
		XEbeta = UtX*Ebeta;

		y_res.setZero();
		y_res = Uty - WEalpha;

		H0 = (sigma2b*D).array()+1;
		H0 = 1/H0.array();

		//////////////////////////////////////
		// full sampling, i.e., sampling effects for all the snps 
		//cout<<tp<<endl;
		size_t tpx=tp;
		if (S%tpx==0)  {
		for (size_t i=0; i<n_snp; i++)   {
			x_col      = UtX.col(i);
			XEbeta    -= x_col*Ebeta(i);

			H = x_col.array()*H0.array();
			tx_ywx_res = H.dot(y_res - XEbeta);

			for (size_t k=0; k<n_k; k++)   {

				if (k == 0) {
					mik_beta(i,k)  = 0; 
					sik2_beta(i,k) = 0; 
					pikexp1(k)     = 0;
					}

				else {
					xtxabk         = H.dot(x_col) + 1/sigma2k(k);
					mik_beta(i,k)  = tx_ywx_res/xtxabk;
					sik2_beta(i,k) = sigma2e/xtxabk;
					pikexp1(k)     = mik_beta(i,k)*mik_beta(i,k)/(2*sik2_beta(i,k)) + log(sqrt(sik2_beta(i,k))) - Elogsigmae - Elogsigmak(k);
				}

				beta_beta(i,k) = gsl_ran_gaussian(rs,sqrt(sik2_beta(i,k))) + mik_beta(i,k);
				pikexp2(k)     = Elogvk(k) + sumElogvl(k);
				}

			index = pikexp1 + pikexp2;
			index = (index.array() - index.maxCoeff()).exp();
			pik_beta.row(i) = index/index.sum();

		// multinomial sampling
		double mult_prob[n_k];
		unsigned int mult_no[n_k];

		for (size_t k=0; k<n_k; k++)   {
			mult_prob[k] = pik_beta(i,k);
			}

		gsl_ran_multinomial(rs, n_k, 1, mult_prob, mult_no);
		for (size_t k=0; k<n_k; k++)   {
			gamma_beta(i,k) = mult_no[k];
		}

		Ebeta(i) = ((beta_beta.row(i)).dot(gamma_beta.row(i)));
		XEbeta  += x_col*Ebeta(i);
		}
		}

		else {
			//partial sampling
			//i.e., sampling effects for pre-selected snps with large effects
			for (size_t t=0; t<s_snp; t++)   {
			  size_t i = snp_no(t);
			x_col      = UtX.col(i);
			XEbeta    -= x_col*Ebeta(i);

			H = x_col.array()*H0.array();
			tx_ywx_res = H.dot(y_res - XEbeta);

			for (size_t k=0; k<n_k; k++)   {

				if (k == 0) {
					mik_beta(i,k)  = 0; 
					sik2_beta(i,k) = 0; 
					pikexp1(k)     = 0;
					}

				else {
					xtxabk         = H.dot(x_col) + 1/sigma2k(k);
					mik_beta(i,k)  = tx_ywx_res/xtxabk;
					sik2_beta(i,k) = sigma2e/xtxabk;
					pikexp1(k)     = mik_beta(i,k)*mik_beta(i,k)/(2*sik2_beta(i,k)) + log(sqrt(sik2_beta(i,k))) - Elogsigmae - Elogsigmak(k);
				}

				beta_beta(i,k) = gsl_ran_gaussian(rs,sqrt(sik2_beta(i,k))) + mik_beta(i,k);
				pikexp2(k)     = Elogvk(k) + sumElogvl(k);
				}

			index = pikexp1 + pikexp2;
			index = (index.array() - index.maxCoeff()).exp();
			pik_beta.row(i) = index/index.sum();

		// multinomial sampling
		double mult_prob[n_k];
		unsigned int mult_no[n_k];

		for (size_t k=0; k<n_k; k++)   {
			mult_prob[k] = pik_beta(i,k);
			}

		gsl_ran_multinomial(rs, n_k, 1, mult_prob, mult_no);
		for (size_t k=0; k<n_k; k++)   {
			gamma_beta(i,k) = mult_no[k];
		}

		Ebeta(i) = ((beta_beta.row(i)).dot(gamma_beta.row(i)));
		XEbeta  += x_col*Ebeta(i);
		}

		}
		//////////////////////////////////////


		if (S > (w_step0-1))   {post_gamma_beta += gamma_beta;}
		//////////////  sampling the mixture snp effects 

		WEalpha.setZero();
		for (size_t j=0;j<n_j; j++)   {
			WEalpha += UtW.col(j)*Ealpha(j);
		}

		y_res.setZero();
		y_res = Uty - XEbeta;
		for (size_t j=0;j<n_j; j++)   {
			WEalpha    -= UtW.col(j)*Ealpha(j);
					  H = UtW.col(j).array()*H0.array();
			m_alpha(j)  = (H.dot(y_res - WEalpha))/(H.dot(UtW.col(j)));
			s2_alpha(j) = sigma2e/(H.dot(UtW.col(j)));
			Ealpha(j)   = m_alpha(j) + gsl_ran_gaussian(rs,sqrt(s2_alpha(j)));
			WEalpha    += UtW.col(j)*Ealpha(j);
			}

		if (S>(w_step0-1))   {post_Ealpha += Ealpha;}

		a_lambda = a0 + n_k;
		b_lambda = b0 - sum_Elogvl2(vk,n_k-1);
		lambdax  = gsl_ran_gamma(rs,a_lambda,1/b_lambda);

		for (size_t k=0; k<n_k-1; k++)   {
			kappa_k(k)  = gamma_beta.col(k).sum() + 1;
			lambda_k(k) = sum_lambda_k(gamma_beta,k,n_k) + lambdax;
			}

		Ebeta2k = beta_beta.cwiseProduct(beta_beta).cwiseProduct(gamma_beta);
		for (size_t k=0; k<n_k; k++)   {
			a_k(k) = (gamma_beta.col(k).sum())/2 + ak;
			b_k(k) = (Ebeta2k.col(k).sum())/(2*sigma2e) + bk(k);
			}

		/*y_res.setZero();
		y_res = Uty - XEbeta - WEalpha;
		for (size_t i=0; i<n_idv; i++) {
			V(i)  = sigma2b*D(i)/(sigma2b*D(i) + 1);
			M(i)  = y_res(i)*V(i);
			Utu(i)= M(i) + gsl_ran_gaussian(rs,sqrt(V(i)*sigma2e));
			if (D(i)  == 0) {
				Ue(i) = 0; 
				Ub(i) = 0;
				}
			else {
				Ue(i) = Utu(i)/(sigma2b*D(i)); //for sigma2e
				Ub(i) = Utu(i)/(sigma2e*D(i)); //for sigma2b
			}
			bv0(i) = y_res(i)*sigma2b/(sigma2b*D(i) + 1); 
			}

		if (S>(w_step0-1))   {bv += bv0;}
		*/

		y_res.setZero();
		y_res = Uty - XEbeta - WEalpha;

		bv0 = (y_res*sigma2b).array()*H0.array();
		if (S>(w_step0-1))   {bv += bv0;}

		H = y_res.array()*H0.array();

		A  = H.dot(y_res);
		B1 = (1/sigma2k.tail(n_k-1).array()).rowwise().replicate(n_snp);
		Ebeta2k   = Ebeta2k.rightCols(n_k-1).cwiseProduct(B1.transpose());
		double B  = Ebeta2k.sum();
		double Gn = gamma_beta.rightCols(n_k-1).sum();
		
		a_e = n_idv/2 + Gn/2 + ae;
		b_e = (A + B + 2*be)/2;

		/*double sigma2bX_new = sigma2bX(S) + gsl_ran_gaussian(rs,sqrt(1));
		double ratio1 = log_sigma2b(D,y_res,sigma2bX_new,sigma2e,ae,be);
		double ratio2 = log_sigma2b(D,y_res,sigma2bX(S), sigma2e,ae,be);
		double ratio  = exp(ratio1-ratio2);
		if (ratio>1) {ratio = 1;}
		else         {ratio = ratio;}
		if (sigma2bX_new<0) {ratio = 0;}
		*/

		//double h_new  = gsl_rng_uniform(rs);
		double h_new  = gsl_ran_beta(rs,2,8);
		double ratio1 = log_h(D,y_res,h_new,sigma2e,ae,be)-log(gsl_ran_beta_pdf(h_new,2,8));
		double ratio2 = log_h(D,y_res,h(S), sigma2e,ae,be)-log(gsl_ran_beta_pdf(h(S),2,8));
		double ratio  = exp(ratio1-ratio2);
		if (ratio>1) {ratio = 1;}
		else         {ratio = ratio;}
		double u = gsl_rng_uniform(rs);
		if (u<ratio) {h(S+1) = h_new;}
		else         {h(S+1) = h(S);}

		sigma2b = h(S+1)/(1-h(S+1));

		//cout<<S<<" "<<setprecision(6)<<h_new<<" "<<sigma2b<<" "<<sigma2e<<" "<<sigma2b*sigma2e<<endl;
		//cout<<setprecision(6)<<sigma2bX_new<<" "<<sigma2b<<" "<<sigma2e<<" "<<sigma2b*sigma2e<<endl;
		if (S > (w_step0-1))   {post_Gn += Gn;}
		////////////////////////////////////

		double llike0=logLike(D,y_res,sigma2b,sigma2e); 
		     llike(S)=llike0; 

		if (S > (w_step0-1))   {	
			post_llike   += llike0;
			post_sigma2e += sigma2e;
			post_sigma2b += sigma2b;
			}

		//size_t xstep = (s_step0 + w_step0)*0.2+1;
		//if (S%xstep==0||S==(s_step0 + w_step0-1)) {
			//ProgressBar ("MCMC sampling ", S, s_step0+w_step0-1);
			//}
		}
	}

	mixture_no = post_Gn/s_step0;
	//cout<<endl<<endl<<"MCMC sampling is finished"<<endl;
	//cout<<endl;

	VectorXd eigen_alpha = VectorXd::Zero (n_snp);
	eigen_alpha = UtX.transpose()*(bv/s_step0)/n_snp;
	//WriteCoeff (eigen_alpha, post_Ebeta/s_step0);
	pheno_mean = post_Ealpha(0)/s_step0;

/////////
	XEbeta.setZero(); //set Gbeta to zero first
		for (size_t i=0; i<n_snp; i++)   {
			x_col   =  UtX.col(i);
			XEbeta += x_col*post_Ebeta(i)/s_step0;
		}
	WEalpha.setZero();
		for (size_t j=0;j<n_j; j++)   {
			WEalpha += UtW.col(j)*post_Ealpha(j)/s_step0;
		}

	y_res.setZero();
	y_res = Uty - XEbeta - WEalpha;

	double llike_hat = logLike(D,y_res,post_sigma2b/s_step0,post_sigma2e/s_step0); 
	VectorXd llike2  = llike.tail(s_step0).array();
	          llike2 = llike2.array() - llike2.sum()/s_step0;
	pD1   =  2*(llike_hat-post_llike/s_step0);
	pD2   =  2*(llike2.dot(llike2))/(s_step0-1);
	if (pD1 < 0)   {pD1 = 1;}
	DIC1  = -2* llike_hat +  2*pD1;
	DIC2  = -2* llike_hat +  2*pD2;
	BIC1  = -2* llike_hat + log(n_idv)*pD1;
	BIC2  = -2* llike_hat + log(n_idv)*pD2;

	//cout<<" K     = " << n_k<<endl;
	//cout<<" pD1   = " << pD1 <<endl;
	//cout<<" pD2   = " << pD2<<endl;
	//cout<<" DIC1  = " << DIC1<<endl;
	//cout<<" DIC2  = " << DIC2<<endl;
	//cout<<" BIC1  = " << BIC1<<endl;
	//cout<<" BIC2  = " << BIC2<<endl;

	//ofstream outfile1;
	//outfile1.open ("llike.txt", ios::out | ios::binary);
	//for (size_t k = 0; k<s_step0; k++)  {
	//	outfile1 << setprecision(7) <<llike2(k)<<endl;
	//	}
	//outfile1 << endl;


	//cout<<endl<< "Computaion Time for DPR.Gibbs = ";
	//cout<<(clock()-time_begin)/(double(CLOCKS_PER_SEC)*60.0)<<" min"<<endl<<endl;
/*
	//output mik and pik
	ofstream outfile1;
	ofstream outfile2;
	ofstream outfile3;
	ofstream outfile4;
	outfile1.open ("mik.txt", ios::out | ios::binary);
	outfile2.open ("pik.txt", ios::out | ios::binary);
	outfile3.open ("sik.txt", ios::out | ios::binary);
	//outfile4.open ("sigmak.txt", ios::out | ios::binary);
	
	for (size_t i = 0; i<n_snp; i++)  {
		for (size_t k = 0; k<n_k; k++)  {
			outfile1 << setprecision(4) <<  mik_beta(i,k)<<"\t";
			outfile2 << setprecision(4) <<  pik_beta(i,k)<<"\t";
			outfile3 << setprecision(4) << sik2_beta(i,k)<<"\t";
			//outfile4 << setprecision(4) << b_k(k)/a_k(k) <<"\t";
			}
		outfile1 << endl;
		outfile2 << endl;
		outfile3 << endl;
		//outfile4 << endl;
	}
	outfile1.close();
	outfile2.close();
	outfile3.close();
	//outfile4.close();
	//output mik and pik
*/

	////////////
	xtx.resize(0);xty.resize(0);wtw.resize(0);wty.resize(0);
	Utu.resize(0);Ue.resize(0);Ub.resize(0);bv0.resize(0);
	bv.resize(0);V.resize(0);M.resize(0);Ealpha.resize(0);
	post_Ealpha.resize(0);m_alpha.resize(0);s2_alpha.resize(0);
	WEalpha.resize(0);XEbeta.resize(0);Ebeta.resize(0);
	post_Ebeta.resize(0);Ebeta2k.resize(0,0);
	B1.resize(0,0);mik_beta.resize(0,0);sik2_beta.resize(0,0);
	pik_beta.resize(0,0);beta_beta.resize(0,0);gamma_beta.resize(0,0);
	sigma2k.resize(0);vk.resize(0);Elogsigmak.resize(0);
	Elogvk.resize(0);sumElogvl.resize(0);pikexp1.resize(0);
	pikexp2.resize(0);index.resize(0);a_k.resize(0);
	b_k.resize(0);kappa_k.resize(0);lambda_k.resize(0);
	y_res.resize(0);bk.resize(0);
	eigen_alpha.resize(0);
}



void LDR::Gibbs_without_u_screen_dic1 (const MatrixXd UtX, const VectorXd Uty, const MatrixXd UtW,const VectorXd D, const VectorXd Wbeta, const VectorXd se_Wbeta, const VectorXd beta, const VectorXi snp_no, double lambda,size_t n_k) 
{
	
	setprecision(5); 
	size_t w_step1=w_step;
	size_t s_step1=s_step;
	clock_t time_begin  = clock();
	size_t n_snp = UtX.cols(); //n * p
	size_t s_snp = snp_no.size();
	size_t n_idv = Uty.size();
	size_t n_j   = UtW.cols();
	VectorXd x_col(n_idv);

	VectorXd xtx(n_snp);
	VectorXd xty(n_snp);
	VectorXd wtw(n_j);
	VectorXd wty(n_j);
	VectorXd Utu = VectorXd::Zero(n_idv);
	VectorXd Ue(n_idv);
	VectorXd Ub(n_idv);

	VectorXd bv0 = VectorXd::Zero(n_idv); //for computation of alpha
	VectorXd bv  = VectorXd::Zero(n_idv); //for computation of alpha
	VectorXd V = VectorXd::Zero (n_idv);//compute Utu 
	VectorXd M = VectorXd::Zero (n_idv);//compute Utu
	VectorXd snp_label(n_snp);

	VectorXd Ealpha      = Wbeta;//intercept
	VectorXd post_Ealpha = VectorXd::Zero(n_j);//save intercept
	VectorXd m_alpha     = Wbeta;
	VectorXd s2_alpha    = se_Wbeta.array()*se_Wbeta.array();

	VectorXd WEalpha(n_idv); 
	double   sigma2e=0, lambdax;
	VectorXd XEbeta(n_idv); //G*beta
	VectorXd Ebeta(n_snp);
	VectorXd post_Ebeta = VectorXd::Zero(n_snp);//save beta
	MatrixXd Ebeta2k;
	MatrixXd B1;
	MatrixXd mik_beta = MatrixXd::Zero(n_snp, n_k);

	VectorXd    llike = VectorXd::Zero(w_step1+s_step1);

	for (size_t i=0;i<n_snp; i++)  {
		//Eigen_matrix_get_row (UtX, i, x_col);
		x_col  = UtX.col(i);
		xtx(i) = x_col.dot(x_col);
		xty(i) = x_col.dot(Uty);
	}

	for (size_t j=0;j<n_j; j++)  {
		wtw(j) = UtW.col(j).dot(UtW.col(j));
		wty(j) = UtW.col(j).dot(Uty);
	}

	//// initial values for SNP and each normal
	//// component has the same initial values
	for (size_t i=0;i<n_snp; i++)  {
		for (size_t k=1;k<n_k; k++)  {
			mik_beta(i,k) = beta(i);
		}
	}

	mik_beta.col(0).fill(0);
	MatrixXd sik2_beta  =(MatrixXd::Random(n_snp, n_k)).array().abs();
	MatrixXd pik_beta   = MatrixXd::Zero  (n_snp, n_k);
	pik_beta = pik_beta.array() + (double)1/n_k;

	MatrixXd beta_beta  = MatrixXd::Zero(n_snp, n_k);
	MatrixXd gamma_beta = MatrixXd::Zero(n_snp, n_k);
	MatrixXd post_gamma_beta = gamma_beta;
	gamma_beta = gamma_beta.array() + (double)1/n_k;

	VectorXd sigma2k (n_k);
	VectorXd vk (n_k);
	VectorXd Elogsigmak (n_k);
	VectorXd Elogvk (n_k);
	VectorXd sumElogvl (n_k);

	VectorXd pikexp1(n_k);
	VectorXd pikexp2(n_k);
	VectorXd index;
	VectorXd a_k = (VectorXd::Random(n_k)).array().abs();
	VectorXd b_k = (VectorXd::Random(n_k)).array().abs();
	VectorXd  kappa_k = pik_beta.colwise().sum();
	VectorXd lambda_k = (VectorXd::Random(n_k)).array().abs();

	XEbeta.setZero();//set to zero
	//for (size_t i=0; i<n_snp; i++)   {
		//x_col   =  UtX.col(i);
		//XEbeta += x_col*beta(i);
		//}
	XEbeta = UtX*Ebeta;

	WEalpha.setZero(); //set to zero
	for (size_t j=0;j<n_j; j++)   {
		WEalpha += UtW.col(j)*m_alpha(j);
		}

	VectorXd y_res  = Uty - WEalpha - XEbeta;
	double sigma2e0 = ((y_res).dot(y_res))/(n_idv-1);
	//cout << endl << endl;
	double a_e  = (2*n_idv + n_snp)/2;
	double b_e  = (a_e-1)*sigma2e0;

	double ak   = 21;
	if (lambda < 0.01)     {lambda =    0.01;}
	else if (lambda > 100) {lambda =     100;}
	else                   {lambda =  lambda;}
	double bk0  = lambda*(ak-1)/n_snp;

	VectorXd bk(n_k);
	bk(0) = bk0;
	for (size_t i=1;i<n_k; i++)  {
		bk(i) = bk(i-1)*1.7*sqrt(pow(i,i));
		}

	for (size_t i=0;i<n_snp; i++)  {
		for (size_t k=0;k<n_k; k++)  {
			sik2_beta(i,k) = bk(k)/(ak-1)*sigma2e0;
		}
	}

	double ae   = 0.1;
	double be   = 0.1;
	//double a0   = 40;
	//double b0   = 4;
	double a0   = 1;
	double b0   = 0.1;
	double a_lambda = a0 + n_k;
	double b_lambda = b0;
	double Elogsigmae, tx_ywx_res, xtxabk, A, post_Gn=0;
	double sigma2b = 0.1;

	double post_llike=0, post_sigma2e=0, post_sigma2b=0;

	VectorXd sigma2bX = (VectorXd::Random(w_step1+s_step1+1)).array().abs();
	sigma2bX(0)=sigma2b;

	VectorXd h = (VectorXd::Random(w_step1+s_step1+1)).array().abs();
	VectorXd H0(n_idv);
	VectorXd H(n_idv);

	//for (size_t k=0; k<(w_step1+s_step1+1); k++)  {
		//cout<<h(k)<<endl;
	//}

	////initial values for a_k, b_k and sigma2k
	Ebeta2k = (mik_beta.cwiseProduct(mik_beta)) + sik2_beta;
	for (size_t k=0; k<n_k; k++)  {
		a_k(k) = (pik_beta.col(k).sum())/2 + ak;
		b_k(k) = (Ebeta2k.col(k).sum())*(a_e/b_e)/2 + bk(k);
		sigma2k(k) = b_k(k)/(a_k(k)-1);
		}

	A = (y_res).dot(y_res);
	B1=(1/sigma2k.array()).rowwise().replicate(n_snp);
	B1.col(0).fill(0);

	Ebeta2k = Ebeta2k.cwiseProduct(B1.transpose());
	double B = Ebeta2k.sum();
	a_e = n_idv + n_snp*0.1 + ae;
	b_e = (A + B + 2*be)/2;

	////random seed
	gsl_rng * rs = gsl_rng_alloc (gsl_rng_mt19937);

	if (n_k == 1)   {
		WritelmmBeta(beta);
		pheno_mean = m_alpha(0);
		}

	// //  begin MCMC sampling
	else {
	for (size_t S=0; S<(w_step1 + s_step1); S++)   {

		sigma2e    = 1/gsl_ran_gamma(rs,a_e,1/b_e);
		Elogsigmae = log(sqrt(sigma2e));

		// save Ebeta for the mixture normal component
		Ebeta = (beta_beta.cwiseProduct(gamma_beta)).rowwise().sum();

		if (S > (w_step1-1))   {post_Ebeta += Ebeta;}

		//sample sigma2k and compute related quantities
		for (size_t k=0; k<n_k; k++)   {
			if (k==0) {
				sigma2k(k)    = 0; 
				Elogsigmak(k) = 0;
				}

			else {
				sigma2k(k)    = 1/gsl_ran_gamma(rs,a_k(k),1/b_k(k));
				Elogsigmak(k) = log(sqrt(sigma2k(k)));
				}

				if (k == (n_k-1)) {vk(k) = 0; Elogvk(k) = 0;}
				else {
					vk(k)     = gsl_ran_beta(rs,kappa_k(k),lambda_k(k));
					Elogvk(k) = log(vk(k));
				}
				sumElogvl(k)  = sum_Elogvl2(vk,k);
			}

		//////////////  sampling the mixture snp effects 
		XEbeta.setZero(); //set Gbeta to zero first
			//for (size_t i=0; i<n_snp; i++)   {
				//x_col   =  UtX.col(i);
				//XEbeta += x_col*Ebeta(i);
		//}
		XEbeta = UtX*Ebeta;

		y_res.setZero();
		y_res = Uty - WEalpha;

		H0 = (sigma2b*D).array()+1;
		H0 = 1/H0.array();

		//////////////////////////////////////
		// full sampling, i.e., sampling effects for all the snps 
		//cout<<tp<<endl;
		size_t tpx=tp;
		if (S%tpx==0)  {
		for (size_t i=0; i<n_snp; i++)   {
			x_col      = UtX.col(i);
			XEbeta    -= x_col*Ebeta(i);

			H = x_col.array()*H0.array();
			tx_ywx_res = H.dot(y_res - XEbeta);

			for (size_t k=0; k<n_k; k++)   {

				if (k == 0) {
					mik_beta(i,k)  = 0; 
					sik2_beta(i,k) = 0; 
					pikexp1(k)     = 0;
					}

				else {
					xtxabk         = H.dot(x_col) + 1/sigma2k(k);
					mik_beta(i,k)  = tx_ywx_res/xtxabk;
					sik2_beta(i,k) = sigma2e/xtxabk;
					pikexp1(k)     = mik_beta(i,k)*mik_beta(i,k)/(2*sik2_beta(i,k)) + log(sqrt(sik2_beta(i,k))) - Elogsigmae - Elogsigmak(k);
				}

				beta_beta(i,k) = gsl_ran_gaussian(rs,sqrt(sik2_beta(i,k))) + mik_beta(i,k);
				pikexp2(k)     = Elogvk(k) + sumElogvl(k);
				}

			index = pikexp1 + pikexp2;
			index = (index.array() - index.maxCoeff()).exp();
			pik_beta.row(i) = index/index.sum();

		// multinomial sampling
		double mult_prob[n_k];
		unsigned int mult_no[n_k];

		for (size_t k=0; k<n_k; k++)   {
			mult_prob[k] = pik_beta(i,k);
			}

		gsl_ran_multinomial(rs, n_k, 1, mult_prob, mult_no);
		for (size_t k=0; k<n_k; k++)   {
			gamma_beta(i,k) = mult_no[k];
		}

		Ebeta(i) = ((beta_beta.row(i)).dot(gamma_beta.row(i)));
		XEbeta  += x_col*Ebeta(i);
		}
		}

		else {
			//partial sampling
			//i.e., sampling effects for pre-selected snps with large effects
			for (size_t t=0; t<s_snp; t++)   {
			  size_t i = snp_no(t);
			x_col      = UtX.col(i);
			XEbeta    -= x_col*Ebeta(i);

			H = x_col.array()*H0.array();
			tx_ywx_res = H.dot(y_res - XEbeta);

			for (size_t k=0; k<n_k; k++)   {

				if (k == 0) {
					mik_beta(i,k)  = 0; 
					sik2_beta(i,k) = 0; 
					pikexp1(k)     = 0;
					}

				else {
					xtxabk         = H.dot(x_col) + 1/sigma2k(k);
					mik_beta(i,k)  = tx_ywx_res/xtxabk;
					sik2_beta(i,k) = sigma2e/xtxabk;
					pikexp1(k)     = mik_beta(i,k)*mik_beta(i,k)/(2*sik2_beta(i,k)) + log(sqrt(sik2_beta(i,k))) - Elogsigmae - Elogsigmak(k);
				}

				beta_beta(i,k) = gsl_ran_gaussian(rs,sqrt(sik2_beta(i,k))) + mik_beta(i,k);
				pikexp2(k)     = Elogvk(k) + sumElogvl(k);
				}

			index = pikexp1 + pikexp2;
			index = (index.array() - index.maxCoeff()).exp();
			pik_beta.row(i) = index/index.sum();

		// multinomial sampling
		double mult_prob[n_k];
		unsigned int mult_no[n_k];

		for (size_t k=0; k<n_k; k++)   {
			mult_prob[k] = pik_beta(i,k);
			}

		gsl_ran_multinomial(rs, n_k, 1, mult_prob, mult_no);
		for (size_t k=0; k<n_k; k++)   {
			gamma_beta(i,k) = mult_no[k];
		}

		Ebeta(i) = ((beta_beta.row(i)).dot(gamma_beta.row(i)));
		XEbeta  += x_col*Ebeta(i);
		}

		}
		//////////////////////////////////////


		if (S > (w_step1-1))   {post_gamma_beta += gamma_beta;}
		//////////////  sampling the mixture snp effects 

		WEalpha.setZero();
		for (size_t j=0;j<n_j; j++)   {
			WEalpha += UtW.col(j)*Ealpha(j);
		}

		y_res.setZero();
		y_res = Uty - XEbeta;
		for (size_t j=0;j<n_j; j++)   {
			WEalpha    -= UtW.col(j)*Ealpha(j);
					  H = UtW.col(j).array()*H0.array();
			m_alpha(j)  = (H.dot(y_res - WEalpha))/(H.dot(UtW.col(j)));
			s2_alpha(j) = sigma2e/(H.dot(UtW.col(j)));
			Ealpha(j)   = m_alpha(j) + gsl_ran_gaussian(rs,sqrt(s2_alpha(j)));
			WEalpha    += UtW.col(j)*Ealpha(j);
			}

		if (S>(w_step1-1))   {post_Ealpha += Ealpha;}

		a_lambda = a0 + n_k;
		b_lambda = b0 - sum_Elogvl2(vk,n_k-1);
		lambdax  = gsl_ran_gamma(rs,a_lambda,1/b_lambda);

		for (size_t k=0; k<n_k-1; k++)   {
			kappa_k(k)  = gamma_beta.col(k).sum() + 1;
			lambda_k(k) = sum_lambda_k(gamma_beta,k,n_k) + lambdax;
			}

		Ebeta2k = beta_beta.cwiseProduct(beta_beta).cwiseProduct(gamma_beta);
		for (size_t k=0; k<n_k; k++)   {
			a_k(k) = (gamma_beta.col(k).sum())/2 + ak;
			b_k(k) = (Ebeta2k.col(k).sum())/(2*sigma2e) + bk(k);
			}

		/*y_res.setZero();
		y_res = Uty - XEbeta - WEalpha;
		for (size_t i=0; i<n_idv; i++) {
			V(i)  = sigma2b*D(i)/(sigma2b*D(i) + 1);
			M(i)  = y_res(i)*V(i);
			Utu(i)= M(i) + gsl_ran_gaussian(rs,sqrt(V(i)*sigma2e));
			if (D(i)  == 0) {
				Ue(i) = 0; 
				Ub(i) = 0;
				}
			else {
				Ue(i) = Utu(i)/(sigma2b*D(i)); //for sigma2e
				Ub(i) = Utu(i)/(sigma2e*D(i)); //for sigma2b
			}
			bv0(i) = y_res(i)*sigma2b/(sigma2b*D(i) + 1); 
			}

		if (S>(w_step1-1))   {bv += bv0;}
		*/

		y_res.setZero();
		y_res = Uty - XEbeta - WEalpha;

		bv0 = (y_res*sigma2b).array()*H0.array();
		if (S>(w_step1-1))   {bv += bv0;}

		H = y_res.array()*H0.array();

		A  = H.dot(y_res);
		B1 = (1/sigma2k.tail(n_k-1).array()).rowwise().replicate(n_snp);
		Ebeta2k   = Ebeta2k.rightCols(n_k-1).cwiseProduct(B1.transpose());
		double B  = Ebeta2k.sum();
		double Gn = gamma_beta.rightCols(n_k-1).sum();
		
		a_e = n_idv/2 + Gn/2 + ae;
		b_e = (A + B + 2*be)/2;

		/*double sigma2bX_new = sigma2bX(S) + gsl_ran_gaussian(rs,sqrt(1));
		double ratio1 = log_sigma2b(D,y_res,sigma2bX_new,sigma2e,ae,be);
		double ratio2 = log_sigma2b(D,y_res,sigma2bX(S), sigma2e,ae,be);
		double ratio  = exp(ratio1-ratio2);
		if (ratio>1) {ratio = 1;}
		else         {ratio = ratio;}
		if (sigma2bX_new<0) {ratio = 0;}
		*/

		//double h_new  = gsl_rng_uniform(rs);
		double h_new  = gsl_ran_beta(rs,2,8);
		double ratio1 = log_h(D,y_res,h_new,sigma2e,ae,be)-log(gsl_ran_beta_pdf(h_new,2,8));
		double ratio2 = log_h(D,y_res,h(S), sigma2e,ae,be)-log(gsl_ran_beta_pdf(h(S),2,8));
		double ratio  = exp(ratio1-ratio2);
		if (ratio>1) {ratio = 1;}
		else         {ratio = ratio;}
		double u = gsl_rng_uniform(rs);
		if (u<ratio) {h(S+1) = h_new;}
		else         {h(S+1) = h(S);}

		sigma2b = h(S+1)/(1-h(S+1));

		//cout<<S<<" "<<setprecision(6)<<h_new<<" "<<sigma2b<<" "<<sigma2e<<" "<<sigma2b*sigma2e<<endl;
		//cout<<setprecision(6)<<sigma2bX_new<<" "<<sigma2b<<" "<<sigma2e<<" "<<sigma2b*sigma2e<<endl;
		if (S > (w_step1-1))   {post_Gn += Gn;}
		////////////////////////////////////

		double llike0=logLike(D,y_res,sigma2b,sigma2e); 
		     llike(S)=llike0; 

		if (S > (w_step1-1))   {	
			post_llike   += llike0;
			post_sigma2e += sigma2e;
			post_sigma2b += sigma2b;
			}

		size_t xstep = (s_step1 + w_step1)*0.2+1;
		if (S%xstep==0||S==(s_step1 + w_step1-1)) {
			ProgressBar ("MCMC sampling ", S, s_step1+w_step1-1);
			}
		}
	}

	mixture_no = post_Gn/s_step1;
	cout<<endl<<"MCMC sampling is finished"<<endl;

	VectorXd eigen_alpha = VectorXd::Zero (n_snp);
	eigen_alpha = UtX.transpose()*(bv/s_step1)/n_snp;
	WriteCoeff (eigen_alpha, post_Ebeta/s_step1);
	pheno_mean = post_Ealpha(0)/s_step1;

/////////
	XEbeta.setZero(); //set Gbeta to zero first
		for (size_t i=0; i<n_snp; i++)   {
			x_col   =  UtX.col(i);
			XEbeta += x_col*post_Ebeta(i)/s_step1;
		}
	WEalpha.setZero();
		for (size_t j=0;j<n_j; j++)   {
			WEalpha += UtW.col(j)*post_Ealpha(j)/s_step1;
		}

	y_res.setZero();
	y_res = Uty - XEbeta - WEalpha;

	double llike_hat = logLike(D,y_res,post_sigma2b/s_step1,post_sigma2e/s_step1); 
	VectorXd llike2  = llike.tail(s_step1).array();
	          llike2 = llike2.array() - llike2.sum()/s_step1;
	pD1   =  2*(llike_hat-post_llike/s_step1);
	pD2   =  2*(llike2.dot(llike2))/(s_step1-1);
	if (pD1 < 0)   {pD1 = 1;}
	DIC1  = -2* llike_hat +  2*pD1;
	DIC2  = -2* llike_hat +  2*pD2;
	BIC1  = -2* llike_hat + log(n_idv)*pD1;
	BIC2  = -2* llike_hat + log(n_idv)*pD2;

	//cout<<" K     = " << n_k<<endl;
	//cout<<" pD1   = " << pD1 <<endl;
	//cout<<" pD2   = " << pD2<<endl;
	//cout<<" DIC1  = " << DIC1<<endl;
	//cout<<" DIC2  = " << DIC2<<endl;
	//cout<<" BIC1  = " << BIC1<<endl;
	//cout<<" BIC2  = " << BIC2<<endl;

	//ofstream outfile1;
	//outfile1.open ("llike.txt", ios::out | ios::binary);
	//for (size_t k = 0; k<s_step1; k++)  {
	//	outfile1 << setprecision(7) <<llike2(k)<<endl;
	//	}
	//outfile1 << endl;


	cout<<"Computaion Time for adaptive DPR.Gibbs = ";
	cout<<(clock()-time_begin)/(double(CLOCKS_PER_SEC)*60.0)<<" min"<<endl<<endl;
/*
	//output mik and pik
	ofstream outfile1;
	ofstream outfile2;
	ofstream outfile3;
	ofstream outfile4;
	outfile1.open ("mik.txt", ios::out | ios::binary);
	outfile2.open ("pik.txt", ios::out | ios::binary);
	outfile3.open ("sik.txt", ios::out | ios::binary);
	//outfile4.open ("sigmak.txt", ios::out | ios::binary);
	
	for (size_t i = 0; i<n_snp; i++)  {
		for (size_t k = 0; k<n_k; k++)  {
			outfile1 << setprecision(4) <<  mik_beta(i,k)<<"\t";
			outfile2 << setprecision(4) <<  pik_beta(i,k)<<"\t";
			outfile3 << setprecision(4) << sik2_beta(i,k)<<"\t";
			//outfile4 << setprecision(4) << b_k(k)/a_k(k) <<"\t";
			}
		outfile1 << endl;
		outfile2 << endl;
		outfile3 << endl;
		//outfile4 << endl;
	}
	outfile1.close();
	outfile2.close();
	outfile3.close();
	//outfile4.close();
	//output mik and pik
*/

	////////////
	xtx.resize(0);xty.resize(0);wtw.resize(0);wty.resize(0);
	Utu.resize(0);Ue.resize(0);Ub.resize(0);bv0.resize(0);
	bv.resize(0);V.resize(0);M.resize(0);Ealpha.resize(0);
	post_Ealpha.resize(0);m_alpha.resize(0);s2_alpha.resize(0);
	WEalpha.resize(0);XEbeta.resize(0);Ebeta.resize(0);
	post_Ebeta.resize(0);Ebeta2k.resize(0,0);
	B1.resize(0,0);mik_beta.resize(0,0);sik2_beta.resize(0,0);
	pik_beta.resize(0,0);beta_beta.resize(0,0);gamma_beta.resize(0,0);
	sigma2k.resize(0);vk.resize(0);Elogsigmak.resize(0);
	Elogvk.resize(0);sumElogvl.resize(0);pikexp1.resize(0);
	pikexp2.resize(0);index.resize(0);a_k.resize(0);
	b_k.resize(0);kappa_k.resize(0);lambda_k.resize(0);
	y_res.resize(0);bk.resize(0);
	eigen_alpha.resize(0);
}

