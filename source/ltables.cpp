#include "ltables.hpp"
#include "grids.hpp"
#include "polyintegration.hpp"

#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <map>
#include <cmath>
#include <complex>
#include <iomanip>

lTables::lTables(int argc, const char *argv[])
{
	m_error = false;

	std::vector<std::string> inputs; for (int i=2; i<argc; i++) inputs.push_back(argv[i]);

	if ((inputs.size() == 1) && (inputs[0] == "-h")) {
		std::cout << "default values: --sNN=5020GeV --pName=Charm --xB=0.6 --LdndxMaxPoints=500000 --LCollMaxPoints=10000 --TCRIT=0.155" << std::endl;
		m_error = true;
	}

	std::map<std::string, std::string> inputparams;
	for (const auto &in : inputs)
	{
 	   	std::string key = in.substr(0, in.find("="));
 	   	std::string::size_type n = 0; while ((n = key.find("-", n)) != std::string::npos) {key.replace(n, 1, ""); n += 0;} //replacing all '-'
		std::string val = in.substr(in.find("=")+1, in.length());
		inputparams[key] = val;
	}

	//checking if configuration file is provided:
	std::map<std::string, std::string> inputparams_f;
	if (inputparams.count("c") > 0) {
		std::ifstream file_in(inputparams["c"]);
		if (!file_in.is_open()) {
			std::cerr << "Error: unable to open configuration file. Aborting..." << std::endl;
			m_error = true;
		}
		std::string line, key, sep, val;
		while (std::getline(file_in, line))
		{
			std::stringstream ss(line);
			ss >> key; ss >>sep; ss >> val;
			inputparams_f[key] = val;
		}
		file_in.close();
	}

	//setting parameter values based on config file values and overwriting with command line values:
	//
	m_sNN = "5020GeV"; if (inputparams_f.count("sNN") > 0) m_sNN = inputparams_f["sNN"];
					   if (inputparams.count("sNN")   > 0) m_sNN =   inputparams["sNN"];

	m_pName = "Charm"; if (inputparams_f.count("pName") > 0) m_pName = inputparams_f["pName"];
					   if (inputparams.count("pName")   > 0) m_pName =   inputparams["pName"];

	m_xB = 0.6; if (inputparams_f.count("xB") > 0) m_xB = stod(inputparams_f["xB"]);
				if (inputparams.count("xB")   > 0) m_xB = stod(  inputparams["xB"]);

	m_LdndxMaxPoints = 500000; if (inputparams_f.count("LdndxMaxPoints") > 0) m_LdndxMaxPoints = stoi(inputparams_f["LdndxMaxPoints"]);
						       if (  inputparams.count("LdndxMaxPoints") > 0) m_LdndxMaxPoints = stoi(  inputparams["LdndxMaxPoints"]);

	m_LCollMaxPoints = 10000; if (inputparams_f.count("LCollMaxPoints") > 0) m_LCollMaxPoints = stoi(inputparams_f["LCollMaxPoints"]);
						      if (inputparams.count("LCollMaxPoints")   > 0) m_LCollMaxPoints = stoi(  inputparams["LCollMaxPoints"]);
	
	m_TCRIT = 0.155; if (inputparams_f.count("TCRIT") > 0) m_TCRIT = stod(inputparams_f["TCRIT"]);
					 if (  inputparams.count("TCRIT") > 0) m_TCRIT = stod(  inputparams["TCRIT"]);

	//checking if provided value of sNN is an option:
	if ((m_sNN != "5440GeV") && (m_sNN != "5020GeV") && (m_sNN != "2760GeV") && (m_sNN != "200GeV")) {
		std::cerr << "Error: provided sNN parameter not an option, please try 5440GeV, 5020GeV, 2760GeV or 200GeV. Aborting..." << std::endl;
		m_error = true;
	}

	m_nf = m_sNN   == "200GeV" ? 2.5 : 3.0;
	m_CR = m_pName == "Gluon"  ? 3.0 : 4.0/3.0;
}

lTables::~lTables() {}

void lTables::runLTables()
{
	if (m_error) return;

	m_Grids.setGridPoints(m_sNN, m_pName, m_TCRIT);

    RadLTables();

    CollLTables();

	if (exportLTables() != 1) return;
}

double lTables::haltonSequence(int index, int base) const
{
	double f = 1.0;
	double res = 0.0;

	while (index > 0) {
		f = f / static_cast<double>(base);
		res += f * static_cast<double>(index % base);
		index = index / base; // integer division
	}

	return res;
}

void lTables::LdndxHSeqInit()
{
	for (size_t i=0; i<m_LdndxMaxPoints; i++) {
		m_LdndxHSeq1.push_back(haltonSequence((i+1)*409, 2));
		m_LdndxHSeq2.push_back(haltonSequence((i+1)*409, 3));
		m_LdndxHSeq3.push_back(haltonSequence((i+1)*409, 5));
	}
}

double lTables::productLog(double x) const
{
	if (x == 0.0) {
		return 0.0;
	}

	double w0, w1;
	if (x > 0.0) {
		w0 = std::log(1.2 * x / std::log(2.4 * x / std::log1p(2.4 * x)));
	}
	else {
		double v = 1.4142135623730950488 * std::sqrt(1.0 + 2.7182818284590452354 * x);
		double N2 = 10.242640687119285146 + 1.9797586132081854940 * v;
		double N1 = 0.29289321881345247560 * (1.4142135623730950488 + N2);
		w0 = -1 + v * (N2 + v) / (N2 + v + N1 * v);
	}

	while (true) {
		double e = std::exp(w0);
		double f = w0 * e - x;
		w1 = w0 - f / ((e * (w0 + 1.0) - (w0 + 2.0) * f / (w0 + w0 + 2.0)));
		if (std::abs(w0 / w1 - 1.0) < 1.4901161193847656e-8) {
			break;
		}
		w0 = w1;
	}
	return w1;
}

double lTables::unitStep(double x) const {
    return (x < 0.0) ? 0.0 : 1.0;
}

long double lTables::unitStep(long double x) const {
    return (x < 0.0L) ? 0.0L : 1.0L;
}

double lTables::dElossDYN(double tau, double p, double x, double k, double q, double varphi, double T) const
{
	double mu = 0.197*std::sqrt((-8.0*(6.0+m_nf)*M_PI*M_PI*T*T)/(2.0*m_nf-33.0)/m_lambda/m_lambda/productLog((-8.0*(6.0 + m_nf)*M_PI*M_PI*T*T)/(2.0*m_nf-33.0)/m_lambda/m_lambda));
	double mg = mu / std::sqrt(2.0);
	double M = 0.0;
	if (m_pName == "Bottom") M = 4.75;
	else if (m_pName == "Charm") M = 1.2;
	else if (m_pName == "Gluon") M = mu/std::sqrt(2.0);
	else M = mu/std::sqrt(6.0);

	double b = std::sqrt(mg*mg + M * M*x*x);
	double e = std::sqrt(p*p + M * M);
	double alpha  = 4.0*M_PI/(11.0 - 2.0*m_nf/3.0)/std::log((k*k + mg*mg + M*M*x*x)/x/m_lambda/m_lambda);
	double alpha1 = 4.0*M_PI/(11.0 - 2.0*m_nf/3.0)/std::log(e*T/0.2/0.2);

	double fn = 1.0;
	fn *= 1.0 / 0.197*m_CR*alpha/M_PI*3.0*alpha1*T*2.0*k*q/M_PI;
	fn *= (mu*mu - mu*mu*m_xB*m_xB)/(q*q + mu*mu*m_xB*m_xB)/(q*q + mu*mu);

	double psi = (k*k + q*q + 2.0*k*q*std::cos(varphi) + b*b)/2.0/x/e*tau/0.197;

	fn *= (1 - std::cos(psi));
	fn *= 2.0/(k*k + b*b)/(k*k + q*q + 2.0*k*q*cos(varphi) + b*b)/(k*k + q*q + 2.0*k*q*std::cos(varphi) + b*b);
	fn *= (-1.0*k*q*std::cos(varphi)*(k*k + q*q + 2.0*k*q*std::cos(varphi)) + b*b*(k*q*std::cos(varphi) + q*q));

	return fn;
}

double lTables::Ldndx(double tau, double p, double T, double x) const
{
	double mu = 0.197*std::sqrt((-8.0*(6.0+m_nf)*M_PI*M_PI*T*T)/(2.0*m_nf-33.0)/m_lambda/m_lambda/productLog((-8.0*(6.0+m_nf)*M_PI*M_PI*T*T)/(2.0*m_nf-33.0)/m_lambda/m_lambda));
	double mg = mu / std::sqrt(2.0);
	double M = 0.0;
	if (m_pName == "Bottom") M = 4.75;
	else if (m_pName == "Charm") M = 1.2;
	else if (m_pName == "Gluon") M = mg;
	else M = mu/sqrt(6.0);
	double e = sqrt(p*p + M * M);

	double kl = 0.00000001; 
	double kh = 2.0*x*(1 - x)*e;
	double kq = (kh - kl);
	double ql = 0.000001;
	double qh = sqrt(4.0*e*T);
	double qq = qh - ql;
	double phil = 0.0;
	double phih = M_PI;
	double phiq = (phih - phil);
	double sum = 0.0; //integration sum
	double k, q, phi; //integration variables

	#pragma omp parallel for reduction(+:sum) private(k,q,phi)
	for (size_t i = 0; i<m_LdndxMaxPoints; i++) {
		  k  =   kl + m_LdndxHSeq1[i]*kq;
		  q  =   ql + m_LdndxHSeq2[i]*qq;
		phi  = phil + m_LdndxHSeq3[i]*phiq;
		sum += 2*dElossDYN(tau, p, x, k, q, phi, T)/x;
	}

	return (sum*kq*qq*phiq/static_cast<double>(m_LdndxMaxPoints));
}

void lTables::RadLTables()
{
	LdndxHSeqInit();

	m_LdndxTbl.resize(m_Grids.tauPtsLength(), std::vector<std::vector<std::vector<double>>>(m_Grids.pPtsLength(), std::vector<std::vector<double>>(m_Grids.TPtsLength(), std::vector<double>(m_Grids.xPtsLength(), 0.0))));;
	m_LNormTbl.resize(m_Grids.tauPtsLength(), std::vector<std::vector<double>>(m_Grids.pPtsLength(), std::vector<double>(m_Grids.TPtsLength(), 0.0)));

	double tau, p, T, x, mu, M, xIntegLimitLow, xIntegLimitHigh;

	for (size_t itau=0; itau<m_Grids.tauPtsLength(); itau++) {
		tau = m_Grids.tauPts(itau);

		for (size_t ip=0; ip<m_Grids.pPtsLength(); ip++) {
			p = m_Grids.pPts(ip);

			for (size_t iT=0; iT<m_Grids.TPtsLength(); iT++) {
				T = m_Grids.TPts(iT);
				
				for (size_t ix=0; ix<m_Grids.xPtsLength(); ix++) {
					x = m_Grids.xPts(ix);
					m_LdndxTbl[itau][ip][iT][ix] = Ldndx(tau, p, T, x);
				}
				
				mu = 0.197*std::sqrt((-8.0*(6.0+m_nf)*M_PI*M_PI*T*T)/(2.0*m_nf-33.0)/m_lambda/m_lambda/productLog((-8.0*(6.0+m_nf)*M_PI*M_PI*T*T)/(2.0*m_nf-33.0)/m_lambda/m_lambda));
				if (m_pName == "Bottom") M = 4.75;
				else if (m_pName == "Charm") M = 1.2;
				else if (m_pName == "Gluon") M = mu/sqrt(2.0);
				else M = mu/sqrt(6.0);
				
				xIntegLimitLow = mu/std::sqrt(2.0)/(p + std::sqrt(p*p + M*M));
				if (m_pName == "Gluon") xIntegLimitHigh = 0.5;
				else xIntegLimitHigh = 1.0 - M/(std::sqrt(p*p + M*M) + p);

				m_LNormTbl[itau][ip][iT] = poly::cubicIntegrate(m_Grids.xPts(), m_LdndxTbl[itau][ip][iT], xIntegLimitLow, xIntegLimitHigh);
			}
		}
	}
}

void lTables::LCollHSeqInit()
{
	for (size_t i=0; i<m_LCollMaxPoints; i++) {
		m_LCollHSeq1.push_back(haltonSequence((i+1)*409, 2));
		m_LCollHSeq2.push_back(haltonSequence((i+1)*409, 3));
		m_LCollHSeq3.push_back(haltonSequence((i+1)*409, 5));
	}
}

std::complex<double> lTables::deltaL2(double q, double w, double T) const
{
	double mu = 0.197*sqrt((-8.0*(6.0+m_nf)*M_PI*M_PI*T*T)/(2.0*m_nf-33.0)/m_lambda/m_lambda/productLog((-8.0*(6.0+m_nf)*M_PI*M_PI*T*T)/(2.0*m_nf-33.0)/m_lambda/m_lambda));

	std::complex<double> q_c = q, w_c = w;
	std::complex<double> log_c = std::log((q_c + w_c)/(q_c - w_c));

	std::complex<double> fn = q*q + mu*mu*(1.0 - w/2.0/q*log_c);
	fn  = fn*fn;
	fn += (M_PI*M_PI*mu*mu*mu*mu/4.0*w*w/q/q);

	return (1.0/fn);
}

std::complex<double> lTables::deltaT2(double q, double w, double T) const
{
	double mu = 0.197*sqrt((-8.0*(6.0+m_nf)*M_PI*M_PI*T*T)/(2.0*m_nf-33.0)/m_lambda/m_lambda/productLog((-8.0*(6.0+m_nf)*M_PI*M_PI*T*T)/(2.0*m_nf-33.0)/m_lambda/m_lambda));

	std::complex<double> q_c = q, w_c = w;
	std::complex<double> log_c = std::log((q_c + w_c)/(q_c - w_c));

	std::complex<double> fn = w*w/q/q + w*(q*q - w*w)/2.0/q/q/q*log_c;
	fn *= (mu*mu/2.0);
	fn += (q*q - w*w);
	fn = fn*fn;
	fn += (M_PI*M_PI*mu*mu*mu*mu/4.0*w*w/q/q*(q*q - w*w)*(q*q - w*w)/4.0/q/q/q/q);

	return (1.0/fn);
}

double lTables::ENumFinite(double p, double T) const
{
	double mu = 0.197*sqrt((-8.0*(6.0+m_nf)*M_PI*M_PI*T*T)/(2.0*m_nf-33.0)/m_lambda/m_lambda/productLog((-8.0*(6.0+m_nf)*M_PI*M_PI*T*T)/(2.0*m_nf-33.0)/m_lambda/m_lambda));
	double M = 1.0;
	if (m_pName == "Bottom") M = 4.75;
	else if (m_pName == "Charm") M = 1.2;
	else if (m_pName == "Gluon") M = mu/std::sqrt(2.0);
	else M = mu/std::sqrt(6.0);
	double e = std::sqrt(p*p + M * M);
	double v = p/e;
	double alpha1 = 4.0*M_PI/(11.0 - 2.0/3.0*m_nf)/std::log(e*T/0.2/0.2);
	double alpha2 = 2.0*M_PI/(11.0 - 2.0/m_nf*3.0)/std::log(mu/0.2);

	//ENumFinite1 integral:
	double ENumFiniteSum1 = 0.0;

	double nfCol;

	double k;
	double kl = 0.0001;
	double kh = m_kmaxColl;
	double kq = kh - kl;

	double ql = 0.0001;
	double qh, qq, q, qmaxCol, qh1, qh2;

	double wl, wh, wq, w;

	#pragma omp parallel for reduction(+:ENumFiniteSum1) private(k,nfCol, qh,qq,q,qmaxCol, wl,wh,wq,w)
	for (size_t i=0; i<m_LCollMaxPoints; i++) {
		std::complex<double> fn_comp;

		k = kl + m_LCollHSeq1[i]*kq;
		nfCol = m_Ng/(std::exp(k/T) - 1.0) + m_nf/(std::exp(k/T) + 1.0);

		qmaxCol = std::sqrt(6.0*e*T);
		qh = ((qmaxCol < k) ? qmaxCol : k);
		qq = qh - ql;
		q = ql + m_LCollHSeq2[i]*qq;

		wl = -q;
		wh = q;
		wq = wh - wl;
		w = wl + m_LCollHSeq3[i]*wq;

		fn_comp  = 2.0/0.197*m_CR*alpha1*alpha2/M_PI/v/v*nfCol*w*unitStep(v*v*q*q - w*w);
		fn_comp *= (deltaL2(q, w, T)*((2.0*k + w)*(2.0*k + w) - q*q)/2.0 + deltaT2(q, w, T)*(q*q - w*w)/4.0/q/q/q/q*((2.0*k + w)*(2.0*k + w) + q*q)*(v*v*q*q - w*w));

		ENumFiniteSum1 += fn_comp.real()*qq*wq;
	}

	ENumFiniteSum1 = ENumFiniteSum1*kq/static_cast<double>(m_LCollMaxPoints);
	
	//ENumFinite2 integral:
	double ENumFiniteSum2 = 0.0;

	#pragma omp parallel for reduction(+:ENumFiniteSum2) private(k,nfCol, ql,qh,qh1,qh2,qq,q,qmaxCol, wl,wh,wq,w)
	for (size_t i=0; i<m_LCollMaxPoints; i++) {
		std::complex<double> fn_comp;

		k = kl + m_LCollHSeq1[i]*kq;
		nfCol = m_Ng/(std::exp(k/T) - 1.0) + m_nf/(std::exp(k/T) + 1.0);

		qmaxCol = std::sqrt(6.0*e*T);
		ql = ((qmaxCol < k) ? qmaxCol : k);
		qh1 = 2.0*k*(1.0 - k/e)/(1.0 - v + 2.0*k/e);
		qh2 = ((k > qh1) ? k : qh1);
		qh = ((qmaxCol < qh2) ? qmaxCol : qh2);
		qq = qh - ql;
		q = ql + m_LCollHSeq2[i]*qq;

		wl = q - 2.0*k;
		wh = q;
		wq = wh - wl;
		w = wl + m_LCollHSeq3[i] * wq;

		fn_comp  = 2.0/0.197*m_CR*alpha1*alpha2/M_PI/v/v*nfCol*w*unitStep(v*v*q*q - w*w);
		fn_comp *= (deltaL2(q, w, T)*((2.0*k + w)*(2.0*k + w) - q*q)/2.0 + deltaT2(q, w, T)*(q*q - w*w)/4.0/q/q/q/q*((2.0*k + w)*(2.0*k + w) + q*q)*(v*v*q*q - w*w));

		ENumFiniteSum2 += fn_comp.real()*qq*wq;
	}
	
	ENumFiniteSum2 = ENumFiniteSum2*kq/static_cast<double>(m_LCollMaxPoints);

	return (ENumFiniteSum1 + ENumFiniteSum2);
}

void lTables::CollLTables()
{
	LCollHSeqInit();

	m_LCollTbl.resize(m_Grids.pCollPtsLength(), std::vector<double>(m_Grids.TCollPtsLength(), 0.0));

	for (size_t ip=0; ip<m_Grids.pCollPtsLength(); ip++) {
		for (size_t iT=0; iT<m_Grids.TCollPtsLength(); iT++) {
			m_LCollTbl[ip][iT] = ENumFinite(m_Grids.pCollPts(ip), m_Grids.TCollPts(iT));
		}
	}
}

int lTables::exportLTables() const
{
	std::stringstream xBss; xBss << std::fixed << std::setprecision(1) << m_xB;
	std::stringstream nfss; nfss << std::fixed << std::setprecision(1) << m_nf;

	{//exporting Ldndx table
		const std::string path_out = "./ltables/ldndx_nf=" + nfss.str() + "_" + m_pName + "_xB=" + xBss.str() + ".dat";
		std::ofstream file_out(path_out, std::ios_base::out);
		if (!file_out.is_open()) {
			std::cerr << "Error: unable to open Ldndx table export file." << std::endl;
			return -1;
		}

		file_out << "#";
		file_out << std::fixed << std::setw(12) <<   "tau" << " ";
		file_out << std::fixed << std::setw(14) <<     "p" << " ";
		file_out << std::fixed << std::setw(12) <<     "T" << " ";
		file_out << std::fixed << std::setw(12) <<     "x" << " ";
		file_out << std::fixed << std::setw(17) << "Ldndx" << "\n";

		for (size_t itau=0; itau<m_Grids.tauPtsLength(); itau++) {
			for (size_t ip=0; ip<m_Grids.pPtsLength(); ip++) {
				for (size_t iT=0; iT<m_Grids.TPtsLength(); iT++) {
					for (size_t ix=0; ix<m_Grids.xPtsLength(); ix++) {
						file_out << std::fixed 		<< std::setw(13) << std::setprecision(10) << m_Grids.tauPts(itau) << " ";
						file_out << std::fixed 		<< std::setw(14) << std::setprecision(10) << m_Grids.pPts(ip) << " ";
						file_out << std::fixed 		<< std::setw(12) << std::setprecision(10) << m_Grids.TPts(iT) << " ";
						file_out << std::fixed 		<< std::setw(12) << std::setprecision(10) << m_Grids.xPts(ix) << " ";
						file_out << std::scientific << std::setw(17) << std::setprecision(10) << m_LdndxTbl[itau][ip][iT][ix] << "\n";
					}
				}
			}
		}

		file_out.close();
	}

	{//exporting LNorm table
		const std::string path_out = "./ltables/lnorm_nf=" + nfss.str() + "_" + m_pName + "_xB=" + xBss.str() + ".dat";
		std::ofstream file_out(path_out, std::ios_base::out);
		if (!file_out.is_open()) {
			std::cerr << "Error: unable to open LNorm table export file." << std::endl;
			return -2;
		}

		file_out << "#";
		file_out << std::fixed << std::setw(12) <<   "tau" << " ";
		file_out << std::fixed << std::setw(14) <<     "p" << " ";
		file_out << std::fixed << std::setw(12) <<     "T" << " ";
		file_out << std::fixed << std::setw(17) << "LNorm" << "\n";

		for (size_t itau=0; itau<m_Grids.tauPtsLength(); itau++) {
			for (size_t ip=0; ip<m_Grids.pPtsLength(); ip++) {
				for (size_t iT=0; iT<m_Grids.TPtsLength(); iT++) {
					file_out << std::fixed 		<< std::setw(13) << std::setprecision(10) << m_Grids.tauPts(itau) << " ";
					file_out << std::fixed 		<< std::setw(14) << std::setprecision(10) << m_Grids.pPts(ip) << " ";
					file_out << std::fixed 		<< std::setw(12) << std::setprecision(10) << m_Grids.TPts(iT) << " ";
					file_out << std::scientific << std::setw(17) << std::setprecision(10) << m_LNormTbl[itau][ip][iT] << "\n";
				}
			}
		}

		file_out.close();
	}

	{//exporting LColl table
		std::string path_out = "./ltables/lcoll_nf=" + nfss.str() + "_" + m_pName + ".dat";
		std::ofstream file_out(path_out, std::ios_base::out);
		if (!file_out.is_open()) {
			std::cerr << "Error: unable to open LColl table export file." << std::endl;
			return -3;
		}

		file_out << "#";
		file_out << std::fixed << std::setw(13) <<     "p" << " ";
		file_out << std::fixed << std::setw(12) <<     "T" << " ";
		file_out << std::fixed << std::setw(17) << "LColl" << "\n";

		for (size_t ip=0; ip<m_Grids.pCollPtsLength(); ip++) {
			for (size_t iT=0; iT<m_Grids.TCollPtsLength(); iT++) {
				file_out << std::fixed 		<< std::setw(14) << std::setprecision(10) << m_Grids.pCollPts(ip) << " ";
				file_out << std::fixed 		<< std::setw(12) << std::setprecision(10) << m_Grids.TCollPts(iT) << " ";
				file_out << std::scientific << std::setw(17) << std::setprecision(10) << m_LCollTbl[ip][iT] << "\n";
			}
		}

		file_out.close();
	}

	return 1;
}