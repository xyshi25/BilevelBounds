#pragma warning(disable : 4996)
#include<fstream>
#include<iostream>
#include<string>
#include<sstream>
#include<iterator>
#include<vector>
#include<random>
#include<ctime>
#include <algorithm>
#include <numeric>
#include<ilcplex/ilocplex.h>

#define MAX std::numeric_limits<int>::max()
//using namespace std;
ILOSTLBEGIN


class BKP
{
public:
	BKP()
		: n_(0), ubudget_(0), lbudget_(0), mu_(0), cost_(env), uweight_(env), lweight_(env) {}
	BKP(string filename, string file_path, string log_path) {
		cost_ = IloIntArray(env);
		uweight_ = IloIntArray(env);
		lweight_ = IloIntArray(env);
		log_path_ = log_path;
		if (!ReadFile(filename, file_path)) {
			return;
		}
	}
	~BKP() {
		Clear();
	}

	void Clear() {
		n_ = 0;
		ubudget_ = 0;
		lbudget_ = 0;

		cost_.clear();
		uweight_.clear();
		lweight_.clear();
		env.end();
	}

	bool ReadFile(string filename, string file_path) {
		filename_ = filename;
		ifstream in_file(file_path + filename + ".txt");
		//	Clear();

		if (in_file.is_open())
		{
			string line;
			int line_num = 0;

			while (getline(in_file, line)) {
				line_num++;
				istringstream split(line);
				if (line_num == 1)
				{
					split >> n_ >> ubudget_ >> lbudget_;
					continue;
				}

				vector<string> split_line((istream_iterator<string>(split)), istream_iterator<string>());

				if (split_line.size() == 0) continue;
				if (split_line[0] == "uweight") {
					for (size_t i = 1; i < split_line.size(); i++) uweight_.add(IloInt(stoi(split_line[i])));
					assert(uweight_.getSize() == n_);
				}
				else if (split_line[0] == "lweight")
				{
					for (size_t i = 1; i < split_line.size(); i++) lweight_.add(stoi(split_line[i]));
					assert(lweight_.getSize() == n_);
				}
				else if (split_line[0] == "cost")
				{
					for (size_t i = 1; i < split_line.size(); i++) cost_.add(stoi(split_line[i]));
					assert(cost_.getSize() == n_);
				}

			}

			ComputeMu();
			return 1;
		}
		else
		{
			cout << "Can't open file " << filename << endl;
			return 0;
		}
	}
	
	bool WriteToMps(string filename) {
		string mps_file = "Dataset/DNegMibs/" + filename + ".mps";
		string aux_file = "Dataset/DNegMibs/" + filename + ".txt";
		ofstream mps_out(mps_file);


		mps_out << "NAME\t" << filename << "\n";
		//ROWS
		mps_out << "ROWS\n N \t OBJROW\n L\t R0 \n";

		// COLUMNS
		mps_out << "COLUMNS\n";

		for (size_t i = 0; i < n_; i++)
		{
			mps_out << ' ' << create_name("C", i) << "\t OBJROW \t" << -cost_[i] << "\t R0 \t" << lweight_[i] << "\t\n";
		}

		// RHS
		mps_out << "RHS\n RHS\t R0 \t " << lbudget_ << "\t\n";

		// BOUNDS
		mps_out << "BOUNDS\n";

		for (size_t i = 0; i < n_; i++)
		{
			mps_out << " BV BOUND\t" << create_name("C", i) << "\t 1" << "\t\n";
		}

		mps_out << "ENDATA\n";
		mps_out.close();
		///////////////////////////////// END mps_out

		ofstream aux_out(aux_file);
		aux_out << "N " << n_ << "\nM " << n_ + 1 << endl;
		for (size_t i = 0; i < n_; i++)
		{
			aux_out << "LC " << n_ + i << endl;
		}

		for (size_t i = 0; i <= n_; i++)
		{
			aux_out << "LR " << 1 + i << endl;
		}

		for (size_t i = 0; i < n_; i++)
		{
			aux_out << "LO " << -cost_[i] << endl;
		}
		aux_out << "OS " << 1 << endl;

		for (size_t i = 0; i < n_; i++)
		{
			aux_out << "IC " << uweight_[i] << endl;
		}

		aux_out << "IB " << ubudget_ << endl;
		aux_out.close();


		return 1;
	}
	
	
	void Solver_CCLW(ofstream &out_result) {
		try
		{
			clock_t start_time = clock();
			// updata p_max, w_max, need more implementation
			vector<double> follower_prefer;
			vector<double> leader_prefer;
			for (size_t i = 0; i < n_; i++)
			{
				follower_prefer.push_back(float(cost_[i]) / float(lweight_[i]));
				leader_prefer.push_back(float(lweight_[i]) / float(uweight_[i]));
			}

			vector<int> idx_follower(n_), idx_leader(n_);

			iota(idx_follower.begin(), idx_follower.end(), 0);
			sort(idx_follower.begin(), idx_follower.end(),
				[&follower_prefer](size_t i1, size_t i2) {return follower_prefer[i1] > follower_prefer[i2]; });

			iota(idx_leader.begin(), idx_leader.end(), 0);
			sort(idx_leader.begin(), idx_leader.end(),
				[&leader_prefer](size_t i1, size_t i2) {return leader_prefer[i1] > leader_prefer[i2]; });



			IloInt sum_lweight(0), sum_uweigt(0);
			for (size_t i = 0; i < n_; i++)
			{

				if (sum_uweigt + uweight_[idx_leader[i]] > ubudget_)
				{
					sum_lweight += IloFloor(IloNum((ubudget_ - sum_uweigt) * lweight_[idx_leader[i]]) / IloNum(uweight_[idx_leader[i]]));
					break;
				}
				sum_lweight += lweight_[idx_leader[i]];
				sum_uweigt += uweight_[idx_leader[i]];
			}


			IloNum cost_max(0), lweight_max(0);
			IloInt sum_lweight2(0);
			for (size_t i = 0; i < n_; i++)
			{
				if ((sum_lweight2 + lweight_[idx_follower[i]] >= lbudget_) && (sum_lweight2 <= lbudget_ + sum_lweight))
				{
					if (cost_max < cost_[idx_follower[i]])	cost_max = cost_[idx_follower[i]];

					if (lweight_max < lweight_[idx_follower[i]])	lweight_max = lweight_[idx_follower[i]];
				}

				if (sum_lweight2 + lweight_[idx_follower[i]] > lbudget_ + sum_lweight)	break;


				sum_lweight2 += lweight_[idx_follower[i]];
			}


			// MIP^1
			IloModel mip_model(env);

			// create variables
			IloNumVarArray xVars(env, n_, 0, 1, IloNumVar::Bool);
			IloNumVarArray yVars(env, n_, 0, 1, IloNumVar::Bool);
			IloNumVarArray uVars(env, n_, 0, IloInfinity);
			IloNumVarArray zVars(env, n_ + 1, 0, IloInfinity);

			// obj
			mip_model.add(IloMinimize(env, IloSum(uVars) + zVars[0] * lbudget_, "mip objective function"));

			// constraints
			mip_model.add(IloScalProd(uweight_, xVars) <= ubudget_);
			for (size_t i = 0; i < n_; i++)
			{
				mip_model.add(zVars[i + 1] - cost_[i] * xVars[i] - uVars[i] <= 0);
				mip_model.add(lweight_[i] * zVars[0] + zVars[i + 1] >= cost_[i]);
			}

			// cplex solver
			string cplex_log = log_path_ + "/cplex_log_" + filename_ + "_CCLW.log";

			IloCplex mip_cplex(mip_model);
			mip_cplex.setParam(IloCplex::TiLim, 600);
			mip_cplex.setParam(IloCplex::Param::Threads, 1);
			ofstream msout_file(cplex_log, std::ios_base::app);
			mip_cplex.setOut(msout_file);

			IloCplex ll_cplex(env);
			ofstream lout_file(cplex_log, std::ios_base::app);
			ll_cplex.setParam(IloCplex::Param::Threads, 1);
			ll_cplex.setOut(lout_file);

			IloNum opt_w = IloInfinity;
			IloNumArray opt_x(env, n_), opt_y(env, n_);


			IloRangeArray ng_cuts(env);
			IloNumArray ng_cuts_bound(env);
			int iter = 0;

			ng_cuts.add(zVars[0] * (lbudget_ - lweight_max) + IloSum(uVars) <= IloSum(cost_));
			ng_cuts_bound.add(0);
			mip_cplex.addCuts(ng_cuts);

			IloInt nodes_max = 0;
			double time_max = 0;
			
			clock_t end_time;
			while (true)
			{
				if ((iter >= 500) && ((int(iter) % 100) == 0))
				{
					end_time = clock();
					if ((float)(end_time - start_time) / CLOCKS_PER_SEC > 3600.0) break;
				}


				IloNumArray iter_x(env, n_), iter_y(env, n_);

				clock_t iter_time = clock();
				if (!mip_cplex.solve())
				{
					//cout << "Infeasible!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
					break;
				}

				//env.out() << "LP relaxation optimal: " << mip_cplex.getObjValue() << endl;

				if (float(clock() - iter_time) / CLOCKS_PER_SEC > time_max) time_max = float(clock() - iter_time) / CLOCKS_PER_SEC;
				
				if (mip_cplex.getNnodes() > nodes_max)	nodes_max = mip_cplex.getNnodes();
				mip_cplex.getValues(iter_x, xVars);

				// stopping critera
				if (opt_w + cost_max <= mip_cplex.getObjValue()) {
					//cout << "OPtimal!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
					break;
				}
				mip_cplex.clearCuts();

				// make maximal of x
				IloNum residule = ubudget_ - IloScalProd(uweight_, iter_x);
				int i = 0;
				while ((i < n_) && (residule > 0))
				{
					if ((iter_x[i] == 0) && (uweight_[i] <= residule))
					{
						residule -= uweight_[i];
						iter_x[i] = 1;
					}
					i++;
				}

				// solve the lower-level problem
				IloModel lower_model(env);
				lower_model.add(yVars);
				lower_model.add(IloScalProd(lweight_, yVars) <= lbudget_);

				// interdiction constraint
				for (int i = 0; i < n_; i++)
				{
					lower_model.add(yVars[i] <= 1 - iter_x[i]);
				}
				//lower_model.add(IloScalProd(iter_x, yVars) == 0);

				// objective for lower level
				lower_model.add(IloMaximize(env, IloScalProd(cost_, yVars), "lower-level objective function"));

				ll_cplex.extract(lower_model);
				//ll_cplex.exportModel("lower_model.lp");

				ll_cplex.solve();
				ll_cplex.getValues(iter_y, yVars);



				if (ll_cplex.getObjValue() < opt_w)
				{
					//cout << "Update Best" << endl;
					opt_w = ll_cplex.getObjValue();
					opt_x = iter_x;
					opt_y = iter_y;

					// update upper bound

					for (size_t i = 0; i < ng_cuts.getSize(); i++)
					{
						ng_cuts[i].setUB(opt_w - ng_cuts_bound[i] - 1);
					}

				}

				//env.out() << iter << '\t' << ll_cplex.getObjValue() << '\t' << iter_x << "\t" << iter_y << endl;
				//mip_cplex.out() << iter << '\t' << opt_w << '\t' << opt_x << "\t" << opt_y << endl;
				//env.out() << iter << '\t' << opt_w << '\t' << opt_x << "\t" << opt_y << endl;
				// add NG cut
				IloExpr new_cut(env);
				IloNum new_cut_bound = ll_cplex.getObjValue();
				for (size_t i = 0; i < n_; i++)
				{
					if (iter_y[i] == 1)
					{
						new_cut += cost_[i] * (-xVars[i]);
					}
				}
				ng_cuts.add(new_cut <= opt_w - new_cut_bound - 1);
				ng_cuts_bound.add(new_cut_bound);

				mip_cplex.addCuts(ng_cuts);

				new_cut.end();
				lower_model.end();
				iter++;

			}

			end_time = clock();
			//mip_cplex.out() << iter << '\t' << opt_w << '\t' << opt_x << "\t" << opt_y << '\t' << (float)(end_time - start_time) / CLOCKS_PER_SEC << endl;
			//env.out() << iter << '\t' << opt_w << '\t' << opt_x << "\t" << opt_y << '\t' << (float)(end_time - start_time) / CLOCKS_PER_SEC << endl;
			out_result << opt_w  << (float)(end_time - start_time) / CLOCKS_PER_SEC << '\t';
			mip_cplex.end();
			ll_cplex.end();
		}
		catch (IloException &e) {
			cerr << "Concert exception caught: " << e << endl;
			return;
		}
		catch (...) {
			cerr << "Unknown exception caught" << endl;
			return;
		}
	}
	bool Solver_local_ext(int K, ofstream &out_result) {
		vector<vector<int> > left_T1, right_T1, left_T2, right_T2;
		vector<int> hweight;
		clock_t start_t = clock();
		// find T for K>=2 
		if (K >= 1)
		{
			local_pair_k(left_T1, right_T1, left_T2, right_T2, hweight, K);
		}
		vector<size_t> idx_h(hweight.size());
		iota(idx_h.begin(), idx_h.end(), 0);

		// sort indexes based on comparing values in v
		sort(idx_h.begin(), idx_h.end(),
			[&hweight](size_t i1, size_t i2) {return hweight[i1] > hweight[i2]; });
		sort(hweight.rbegin(), hweight.rend());

		// reduce the weight size
		vector<int> rhweight; // remove the duplication items in hweight
		vector<size_t> idx_rh(hweight.size()); // the index that connect hweight to rhweight
		for (size_t i = 0; i < hweight.size(); i++) {
			if (i == 0) {
				rhweight.push_back(hweight[i]);
				idx_rh[0] = 0;
				continue;
			}

			if (hweight[i] != hweight[i - 1])	rhweight.push_back(hweight[i]);

			idx_rh[i] = rhweight.size() - 1;
		}

		IloIntArray ilohweight(env, rhweight.size());
		for (size_t i = 0; i < rhweight.size(); i++)
		{
			ilohweight[i] = (IloInt(rhweight[i]));
		}

		try
		{
			IloModel model(env);

			//  create varaibles
			IloNumVarArray xVars(env);
			IloNumVarArray yVars(env);
			IloNumVarArray wVars(env);
			IloNumVarArray vVars(env);
			for (size_t i = 0; i < n_; i++)
			{
				xVars.add(IloNumVar(env, 0, 1, IloNumVar::Bool, create_name("x", i).c_str()));
				yVars.add(IloNumVar(env, 0, 1, IloNumVar::Bool, create_name("y", i).c_str()));
			}
			for (size_t i = 0; i < rhweight.size(); i++)
			{
				wVars.add(IloNumVar(env, 0, 1, IloNumVar::Bool, create_name("w", i).c_str()));
				vVars.add(IloNumVar(env, 0, 1, IloNumVar::Bool, create_name("v", i).c_str()));

			}
			

			model.add(xVars);
			model.add(yVars);
			model.add(wVars);
			model.add(vVars);

			// objective function
			model.add(IloMinimize(env, IloScalProd(cost_, yVars), "objective function"));

			// S_0 constraints
			for (size_t i = 0; i < n_; i++)
			{
				model.add(xVars[i] + yVars[i] <= 1);
			}
			model.add(IloScalProd(uweight_, xVars) <= ubudget_);
			model.add(IloScalProd(lweight_, yVars) <= lbudget_);

			// S_k constraints
			for (size_t i = 0; i < left_T1.size(); i++)
			{
				vector<int> T_1 = left_T1[i], T_2 = right_T1[i];

				IloExpr cons(env);
				for (size_t j = 0; j < T_1.size(); j++) cons += xVars[T_1[j]] + yVars[T_1[j]];
				for (size_t j = 0; j < T_2.size(); j++) cons += 1 - yVars[T_2[j]];

				model.add(cons >= 1);
				cons.end();
			}

			IloRangeArray mixing_set(env);
			IloExpr mixing_con(env);
			mixing_con = IloScalProd(lweight_, yVars);
			for (size_t i = 0; i < rhweight.size(); i++)
			{
				if (i == rhweight.size() - 1)
				{
					mixing_con += (rhweight[i] - mu_) * vVars[i];
				}
				else
				{
					mixing_con += (rhweight[i] - rhweight[i + 1]) * vVars[i];
					mixing_set.add(vVars[i] - vVars[i + 1] >= 0);
				}
				mixing_set.add(wVars[i] - vVars[i] >= 0);
			}
			mixing_set.add(mixing_con >= rhweight[0]);
			model.add(mixing_set);
			mixing_set.end();

			for (size_t i = 0; i < hweight.size(); i++)
			{
				vector<int> T_1 = left_T2[idx_h[i]], T_2 = right_T2[idx_h[i]];
				

				IloExpr con1(env), con2(env);
				for (size_t j = 0; j < T_1.size(); j++) {
					con1 += xVars[T_1[j]] + yVars[T_1[j]];
					con2 += xVars[T_1[j]];
				}
				for (size_t j = 0; j < T_2.size(); j++) {
					con1 += 1 - yVars[T_2[j]];
					con2 += xVars[T_2[j]];
				}

				model.add(con1 - wVars[idx_rh[i]] >= 0);
				//model.add(int(T_1.size() + T_2.size()) * wVars[i] - con2 >= 0);
				con1.end();
				con2.end();
			}

			// cplex solver
			IloCplex cplex(model);
			cplex.setParam(IloCplex::TiLim, 3600);
			cplex.setParam(IloCplex::Param::Threads, 1);
			string cplex_log = log_path_ + "/cplex_log_" + filename_ + "_" + to_string(K) + "_ext.log";
			ofstream out_file(cplex_log, std::ios_base::out);
			cplex.setOut(out_file);
			
			cplex.solve();
			clock_t end_t = clock();;

			IloNum lb = cplex.getObjValue();
			IloNumArray opt_x(env), opt_y(env);
			cplex.getValues(opt_x, xVars);
			//cplex.getValues(opt_y, yVars);

			// solve the lower level problem
			IloModel lower_model(env);
			lower_model.add(yVars);
			lower_model.add(IloScalProd(lweight_, yVars) <= lbudget_);

			// interdiction constraint
			lower_model.add(IloScalProd(opt_x, yVars) == 0);

			// objective for lower level
			lower_model.add(IloMaximize(env, IloScalProd(cost_, yVars), "lower-level objective function"));

			cplex.extract(lower_model);
			cplex.solve();

			IloNum ub = cplex.getObjValue();
			
			double gap = (ub - lb) / float(ub);
			
			cout << K << '\t' << lb << '\t' << ub << '\t' << float(end_t - start_t) / CLOCKS_PER_SEC << endl;
			out_result  << float(end_t - start_t) / CLOCKS_PER_SEC << '\t';

			cplex.end();
			model.end();

			if (lb == ub)
			{
				return 1;
			}
			else
			{
				return 0;
			}

		}
		catch (IloException &e) {
			cerr << "Concert exception caught: " << e << endl;
			return 0;
		}
		catch (...) {
			cerr << "Unknown exception caught" << endl;
			return 0;
		}
	}
	
	// Solver locak with neighborhood size K
	bool Solver_local(int K, ofstream &out_result) {
		vector<vector<int> > left_T1, right_T1, left_T2, right_T2;
		vector<int> hweight;

		clock_t start_t = clock();

		// find T for K>=2 
		if (K >= 1)
		{
			local_pair_k(left_T1, right_T1, left_T2, right_T2, hweight, K);
		}

		vector<size_t> idx_h(hweight.size());
		iota(idx_h.begin(), idx_h.end(), 0);

		// sort indexes based on comparing values in v
		sort(idx_h.begin(), idx_h.end(),
			[&hweight](size_t i1, size_t i2) {return hweight[i1] > hweight[i2]; });
		sort(hweight.rbegin(), hweight.rend());


		// reduce the weight size
		vector<int> rhweight; // remove the duplication items in hweight
		vector<size_t> idx_rh(hweight.size()); // the index that connect hweight to rhweight
		for (size_t i = 0; i < hweight.size(); i++) {
			if (i == 0) {
				rhweight.push_back(hweight[i]);
				idx_rh[0] = 0;
				continue;
			}

			if (hweight[i] != hweight[i - 1])	rhweight.push_back(hweight[i]);

			idx_rh[i] = rhweight.size() - 1;
		}

		IloIntArray ilohweight(env, rhweight.size());
		for (size_t i = 0; i < rhweight.size(); i++)
		{
			ilohweight[i] = (IloInt(rhweight[i]));
		}

		try
		{
			IloModel model(env);

			//  create varaibles
			IloNumVarArray xVars(env);
			IloNumVarArray yVars(env);
			IloNumVarArray wVars(env);
			for (size_t i = 0; i < n_; i++)
			{
				xVars.add(IloNumVar(env, 0, 1, IloNumVar::Bool, create_name("x", i).c_str()));
				yVars.add(IloNumVar(env, 0, 1, IloNumVar::Bool, create_name("y", i).c_str()));
			}
			for (size_t i = 0; i < rhweight.size(); i++)	wVars.add(IloNumVar(env, 0, 1, IloNumVar::Bool, create_name("w", i).c_str()));

			model.add(xVars);
			model.add(yVars);
			model.add(wVars);

			// objective function
			model.add(IloMinimize(env, IloScalProd(cost_, yVars), "objective function"));

			// S_0 constraints
			for (size_t i = 0; i < n_; i++)
			{
				model.add(xVars[i] + yVars[i] <= 1);
			}
			model.add(IloScalProd(uweight_, xVars) <= ubudget_);
			model.add(IloScalProd(lweight_, yVars) <= lbudget_);

			// S_k constraints
			for (size_t i = 0; i < left_T1.size(); i++)
			{
				vector<int> T_1 = left_T1[i], T_2 = right_T1[i];

				IloExpr cons(env);
				for (size_t j = 0; j < T_1.size(); j++) cons += xVars[T_1[j]] + yVars[T_1[j]];
				for (size_t j = 0; j < T_2.size(); j++) cons += 1 - yVars[T_2[j]];

				model.add(cons >= 1);
				cons.end();
			}

			IloRangeArray mixing_set(env);
			for (size_t i = 0; i < rhweight.size(); i++) {
				mixing_set.add(IloScalProd(lweight_, yVars) + (rhweight[i] - mu_) * wVars[i] >= rhweight[i]);
			}
			model.add(mixing_set);
			mixing_set.end();

			for (size_t i = 0; i < hweight.size(); i++)
			{
				vector<int> T_1 = left_T2[idx_h[i]], T_2 = right_T2[idx_h[i]];


				IloExpr con1(env), con2(env);
				for (size_t j = 0; j < T_1.size(); j++) {
					con1 += xVars[T_1[j]] + yVars[T_1[j]];
					con2 += xVars[T_1[j]];
				}
				for (size_t j = 0; j < T_2.size(); j++) {
					con1 += 1 - yVars[T_2[j]];
					con2 += xVars[T_2[j]];
				}

				model.add(con1 - wVars[idx_rh[i]] >= 0);
				con1.end();
				con2.end();
			}



			// cplex solver
			IloCplex cplex(model);
			
			cplex.setParam(IloCplex::TiLim, 3600);
			cplex.setParam(IloCplex::Param::Threads, 1);
			string cplex_log = log_path_ + "/cplex_log_" + filename_ + "_" + to_string(K) + ".log";
			ofstream out_file(cplex_log, std::ios_base::out);
			cplex.setOut(out_file);
	
			cplex.solve();
			clock_t end_t = clock();

			IloNum lb = cplex.getObjValue();
			IloNumArray opt_x(env), opt_y(env), opt_w(env);
			cplex.getValues(opt_x, xVars);
			//cplex.getValues(opt_y, yVars);
			//cplex.getValues(opt_w, wVars);
			//env.out() << "Optimal x: " << opt_x << endl;
			//env.out() << "Optimal y: " << opt_y << endl;
			//env.out() << "Optimal w: " << opt_w << endl;

			// solve the lower level problem
			IloModel lower_model(env);
			lower_model.add(yVars);
			lower_model.add(IloScalProd(lweight_, yVars) <= lbudget_);

			// interdiction constraint
			lower_model.add(IloScalProd(opt_x, yVars) == 0);

			// objective for lower level
			lower_model.add(IloMaximize(env, IloScalProd(cost_, yVars), "lower-level objective function"));

			cplex.extract(lower_model);
			cplex.solve();

			IloNum ub = cplex.getObjValue();

			out_result << left_T1.size() + left_T2.size() << '\t' << rhweight.size() << '\t' << lb << '\t' << ub << '\t' << float(end_t - start_t) / CLOCKS_PER_SEC << '\t';
			cout << left_T1.size() + left_T2.size() << '\t' << rhweight.size() << '\t' << lb << '\t' << ub << '\t' << float(end_t - start_t) / CLOCKS_PER_SEC << endl;

			cplex.end();
			model.end();

			if (lb == ub)
			{
				return 1;
			}
			else
			{
				return 0;
			}

		}
		catch (IloException &e) {
			cerr << "Concert exception caught: " << e << endl;
			return 0;
		}
		catch (...) {
			cerr << "Unknown exception caught" << endl;
			return 0;
		}
	}
	



private:
	void ComputeMu() {

		// compute mu_2
		int mu_2 = lbudget_ - lweight_[n_ - 1];
		
		// compute mu_1
		vector<double> unit_price(n_);
		for (size_t i = 0; i < n_; i++)	unit_price[i] = double(lweight_[i]) / double(uweight_[i]);

		// sort unit price as decreasing order
		vector<int> idx(n_);
		iota(idx.begin(), idx.end(), 0);
		sort(idx.begin(), idx.end(),
			[&unit_price](size_t i1, size_t i2) {return unit_price[i1] > unit_price[i2]; });

		// solve the linear program
		double mu_1 = 0, sum_b = 0;
		for (size_t i = 0; i < n_; i++)
		{
			if (sum_b + uweight_[idx[i]] <= ubudget_)
			{
				mu_1 += lweight_[idx[i]];
				sum_b += uweight_[idx[i]];
			}
			else
			{
				mu_1 += double(ubudget_ - sum_b) * unit_price[idx[i]];
				break;
			}
		}

		double sum_weight = 0;
		for (size_t i = 0; i < n_; i++)	sum_weight += lweight_[i];
		
		mu_1 = ceil(sum_weight - mu_1);

		// finally compute mu_
		mu_ = min(int(mu_1), mu_2);

	}

	string create_name(string letters, int number1, int number2 = -1, int number3 = -1) {
		string name;
		stringstream namestream;
		namestream << letters << number1;
		if (number2 >= 0) {
			namestream << "," << number2;
		}
		if (number3 >= 0) {
			namestream << "," << number3;
		}
		name = namestream.str();
		return name;
	}

	// for BP_k problem, return all pairs whose neighbors are less than K
	void local_pair_k(vector<vector<int> > &left_T1, vector<vector<int> > &right_T1, vector<vector<int> > &left_T2, vector<vector<int> > &right_T2, vector<int> &hweight, int K) {
		vector<int> N;

		// K = 1
		for (size_t i = 0; i < n_; i++)
		{
			N.push_back(i);

			if (lbudget_ + 1 - lweight_[i] < 0)	continue;

			left_T2.push_back(vector<int>(1, i));
			right_T2.push_back(vector<int>());
			hweight.push_back(lbudget_ + 1 - lweight_[i]);
		}

		// K = 2, \ldots, K
		for (size_t k = 2; k <= K; k++)
		{

			for (size_t t = 1; t <= floor(float(k) / 2.0); t++)
			{
				vector<vector<int> > comb = subset_k(N, t);
				for (size_t i = 0; i < comb.size(); i++)
				{
					// |T_1| = t, |T_2| = k-t
					vector<int> T_1 = comb[i];
					vector<int> N_T; // the set complement of T1
					int wsum_T1(0), sum_T1(0), min_T1 = std::numeric_limits<int>::max();
					for (size_t l = 0; l < t; l++)
					{
						wsum_T1 += lweight_[T_1[l]];
						sum_T1 += cost_[T_1[l]];
						if (min_T1 > cost_[T_1[l]]) min_T1 = cost_[T_1[l]];
					}


					set_difference(N.begin(), N.end(), T_1.begin(), T_1.end(),
						std::inserter(N_T, N_T.begin()));

					vector<vector<int> > comb_2 = subset_k(N_T, k - t);
					for (size_t j = 0; j < comb_2.size(); j++)
					{
						vector<int> T_2 = comb_2[j];
						int wsum_T2(0), sum_T2(0), min_T2 = std::numeric_limits<int>::max();
						for (size_t l = 0; l < T_2.size(); l++)
						{
							wsum_T2 += lweight_[T_2[l]];
							sum_T2 += cost_[T_2[l]];
							if (min_T2 > cost_[T_2[l]]) min_T2 = cost_[T_2[l]];
						}

						// T_1 > T_2

						if (lbudget_ + 1 - abs(wsum_T2 - wsum_T1) < 0)	continue;

						if ((sum_T1 > sum_T2) && (sum_T1 - min_T1 <= sum_T2))
						{
							if (wsum_T1 > wsum_T2)
							{
								left_T2.push_back(T_1);
								right_T2.push_back(T_2);
								hweight.push_back(lbudget_ + 1 + wsum_T2 - wsum_T1);
							}
							else
							{
								left_T1.push_back(T_1);
								right_T1.push_back(T_2);
							}

						}
						if ((float(t) != float(k) / 2.0) && (sum_T1 < sum_T2) && (sum_T2 - min_T2 <= sum_T1))
						{
							if (wsum_T1 < wsum_T2)
							{
								left_T2.push_back(T_2);
								right_T2.push_back(T_1);
								hweight.push_back(lbudget_ + 1 - wsum_T2 + wsum_T1);
							}
							else
							{
								left_T1.push_back(T_2);
								right_T1.push_back(T_1);
							}

						}
					}
				}
			}

		}
	}


	// find all subsets of S with size K 
	vector<vector<int> > subset_k(const vector<int> &S, int K) {
		vector<vector<int> > comb;

		if (K == 1)
		{
			for (size_t i = 0; i < S.size(); i++) {
				comb.push_back(vector<int>(1, S[i]));
			}

			return comb;

		}

		vector<int> tmp_S = S;
		for (size_t i = 0; i < S.size(); i++)
		{
			tmp_S.erase(tmp_S.begin());
			vector<vector<int> > tmp_comb = subset_k(tmp_S, K - 1);

			for (size_t j = 0; j < tmp_comb.size(); j++) tmp_comb[j].insert(tmp_comb[j].begin(), S[i]);
			comb.insert(comb.end(), tmp_comb.begin(), tmp_comb.end());

		}

		return comb;
	}

	string filename_;
	string log_path_;
	int n_;
	int ubudget_;
	int lbudget_;
	int mu_;


	IloIntArray cost_;
	IloIntArray uweight_;
	IloIntArray lweight_;
	IloEnv env;

};


int main() {


	/*********************************************************************************************
	*
	*									DNEG
	*
	**********************************************************************************************/

	/*string file_name="DNeg_n20k4j1";
	string file_path = "/home/xus6/BFKO/Dataset/DNeg/";

	ofstream out_file("results.txt");

	BKP test_prob(file_name, file_path);
	test_prob.Solver_local(8, out_file);
	test_prob.Solver_local_ext(1, out_file);
	test_prob.Solver_local_ext(3, out_file);*/

	// Second Experiments
	/*ofstream out_file("DNeg_results.txt");
	string file_path = "/home/xus6/BFKO/Dataset/DNeg/";

	for (int i = 0; i < 10; i++)
	{
		string file_name = "DNeg_n15k4j" + to_string(i);
		cout << file_name << endl;
		BKP test_prob(file_name, file_path);
		out_file << file_name << "\t";
		test_prob.Solver_CCLW(out_file);

		for (int k = 0; k < 11; k++)
		{
			test_prob.Solver_local(k, out_file);
		}
   
    out_file << '\n';

	}*/


	
	ofstream out_file("DNeg_results.txt");
	string file_path = "Dataset/DNeg/";
	string log_path = "Log/Log_KIP";

	out_file << "\t | CCLW | \t \t \t KIP_0 \t \t | \t \t \t KIP_1 \t \t \t |  \t \t \t KIP_2 \t \t \t | \t \t \t KIP_3 \t \t \t |" << endl;
	out_file << "File \t OptObj |\t Time |\t |T| \t #mixing_set \t ObjL \t ObjU \t Time | \t |T| \t #mixing_set \t ObjL \t ObjU \t Time \t ExtTime | \t |T| \t #mixing_set \t ObjL \t ObjU \t Time \t ExtTime | \t |T| \t #mixing_set \t ObjL \t ObjU \t Time \t ExtTime |" << endl;

	for (int n = 10; n <= 50; n+=10)
	{
		for (int k = 0; k < 10; k++)
		{
			for (int i = 0; i < 10; i++)
			{
				string file_name = "DNeg_n" + to_string(n) + "k" + to_string(k) + "j" + to_string(i);

				cout << file_name << endl;
				BKP test_prob(file_name, file_path, log_path);

				out_file << file_name << "\t";
				test_prob.Solver_CCLW(out_file);
				test_prob.Solver_local(0, out_file);
				test_prob.Solver_local(1, out_file);
				test_prob.Solver_local_ext(1, out_file);
				test_prob.Solver_local(2, out_file);
				test_prob.Solver_local_ext(2, out_file);
				test_prob.Solver_local(3, out_file);
				test_prob.Solver_local_ext(3, out_file);
				out_file << '\n';
			}
		}
	}
	
	out_file.close();

	

	// ifstream in_file("test.txt");
	// if (in_file.is_open())
	// {
	// 	string file_name;
	// 	while (getline(in_file, file_name))
	// 	{
	// 		out_file << file_name << "\t";
	// 		cout << file_name << endl;
	// 		BKP test_prob(file_name, file_path);

	// 		test_prob.Solver_CCLW(out_file);
	// 		test_prob.Solver_local(0, out_file);
	// 		test_prob.Solver_local(1, out_file);
	// 		test_prob.Solver_local_ext(1, out_file);
	// 		test_prob.Solver_local(2, out_file);
	// 		test_prob.Solver_local_ext(2, out_file);
	// 		test_prob.Solver_local(3, out_file);
	// 		test_prob.Solver_local_ext(3, out_file);
	// 		out_file << '\n';
	// 	}

	// 	in_file.close();
	// 	out_file.close();
	// }
	return 0;
}