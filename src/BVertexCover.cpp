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
#include <chrono>

#define MAX std::numeric_limits<int>::max()
//using namespace std;
ILOSTLBEGIN

// 3D variable array
typedef IloArray<IloNumVarArray> NumVarMatrix;



class BVC
{
public:
	BVC()
		: budget_(0), nnode_(0), nedge_(0), leader_weight_(env), follower_weight_(env) {}
	BVC(string filename, string file_path, string log_path, bool flag) {
		log_path_ = log_path;
		leader_weight_ = IloIntArray(env);
		follower_weight_ = IloIntArray(env);
		if (!ReadFile(filename, file_path, flag)) {
			cout << "Fail to read file " << filename << endl;
			return;
		}
	}
	~BVC() {
		Clear();
	}

	void Clear() {
		budget_ = 0;
		nnode_ = (0);
		nedge_ = (0);

		leader_weight_.clear();
		follower_weight_.clear();

		env.end();
	}

	bool ReadFile(string filename, string file_path, bool flag) {
		leader_weight_ = IloIntArray(env);
		follower_weight_ = IloIntArray(env);

		filename_ = filename;

		ifstream in_file(file_path + filename + ".txt");

		if (in_file.is_open())
		{
			string line;
			int line_num = 0;

			while (getline(in_file, line))
			{
				line_num++;
				istringstream split(line);
				if (line_num == 1)
				{
					split >> nnode_ >> nedge_ >> budget_;
					neighborhood_.resize(nnode_);
					for (size_t i = 0; i < nnode_; i++)	neighborhood_[i].push_back(i); // initialize neighborhood vector with itself first

					continue;
				}

				vector<string> split_line((istream_iterator<string>(split)), istream_iterator<string>());

				if (split_line.size() == 0) continue;
				if (split_line[0] == "uweight")
				{
					for (size_t i = 1; i < split_line.size(); i++) leader_weight_.add(IloInt(stoi(split_line[i])));
				}
				else if (split_line[0] == "lweight")
				{
					for (size_t i = 1; i < split_line.size(); i++) follower_weight_.add(IloInt(stoi(split_line[i])));
				}
				else if (split_line[0] == "edge")
				{
					for (int edge = 0; edge < nedge_; edge++)
					{
						string edge_line;
						getline(in_file, edge_line);
						istringstream edge_split(edge_line);

						vector<int> vedge(2);
						edge_split >> vedge[0] >> vedge[1];
						vedge_.push_back(vedge);
						neighborhood_[vedge[0]].push_back(vedge[1]);
						neighborhood_[vedge[1]].push_back(vedge[0]);
					}

				}
			}


			// sort neighborhood
			for (int i = 0; i < nnode_; i++)
			{
				sort(neighborhood_[i].begin(), neighborhood_[i].end());
			}

			// interdiction setting
			if (flag)
			{
				for (int i = 0; i < nnode_; i++)
				{
					leader_weight_[i] = follower_weight_[i];
				}

			}

			return 1;
		}
		else
		{
			return 0;
		}
	}

	bool Solver_local(int K, ofstream &out_result) {
		/**********************************************************
		*
		*				Find T^k
		*
		***********************************************************/
		vector<vector<int> > left_T1, right_T2;
		vector<vector<int> > hweight;

		vector<vector<size_t> > index_h, index_t;
		vector<vector<int> > rhweight(nnode_), index_rh;

		clock_t start_t = clock();
		std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

		// find T for K>=2 
		if (K >= 1)
		{
			local_pair_k(left_T1, right_T2, hweight, K);

			for (size_t i = 0; i < nnode_; i++)
			{
				vector<int> hweight_i = hweight[i];
				vector<size_t> idx_i(hweight[i].size()), idx_ti(hweight[i].size());
				iota(idx_i.begin(), idx_i.end(), 0);
				iota(idx_ti.begin(), idx_ti.end(), 0);

				// sort indexes based on comparing values in v
				sort(idx_i.begin(), idx_i.end(),
					[&hweight_i](size_t i1, size_t i2) {return hweight_i[i1] < hweight_i[i2]; });
				sort(idx_ti.begin(), idx_ti.end(),
					[&idx_i](size_t i1, size_t i2) {return idx_i[i1] < idx_i[i2]; });

				sort(hweight[i].begin(), hweight[i].end());
				index_h.push_back(idx_i);
				index_t.push_back(idx_ti);

				// reduce the weight size
				vector<int> rhweight_i; // remove the duplication items in hweight
				vector<int> idx_ri(hweight[i].size(), -1); // the index that connect hweight to rhweight
				for (size_t t = 0; t < hweight[i].size(); t++) {
					if (hweight[i][t] <= 0) continue;

					if (t == 0) {
						rhweight_i.push_back(hweight[i][t]);
						idx_ri[0] = 0;
						continue;
					}

					if (hweight[i][t] != hweight[i][t - 1])	rhweight_i.push_back(hweight[i][t]);

					idx_ri[t] = rhweight_i.size() - 1;
				}

				index_rh.push_back(idx_ri);
				rhweight[i] = (rhweight_i);
			}
		}

		/**********************************************************
		*
		*				Modeling
		*
		***********************************************************/


		try
		{
			IloModel model(env);

			// create varaibles
			IloNumVarArray xVars(env);
			IloNumVarArray yVars(env);
			NumVarMatrix zVars(env, nnode_);
			for (size_t i = 0; i < nnode_; i++)
			{
				xVars.add(IloNumVar(env, 0, 1, IloNumVar::Bool, create_name("x", i).c_str()));
				yVars.add(IloNumVar(env, 0, 1, IloNumVar::Bool, create_name("y", i).c_str()));

				zVars[i] = IloNumVarArray(env);
				for (size_t t = 0; t < rhweight[i].size(); t++)
				{
					//zVars[i].add(IloNumVar(env, 0, 1, IloNumVar::Bool, create_name("z", i, t).c_str()));
					zVars[i].add(IloNumVar(env, 0, 1, IloNumVar::Bool, create_name("z", i, t).c_str()));
				}
				model.add(zVars[i]);
			}

			model.add(xVars);
			model.add(yVars);

			// objective function
			model.add(IloMaximize(env, IloScalProd(leader_weight_, yVars), "Objective"));

			// S_0 constraints
			// leader budget
			model.add(IloSum(xVars) <= budget_);
			IloRangeArray neighbor_cons(env);
			for (size_t i = 0; i < nnode_; i++)
			{
				model.add(xVars[i] + yVars[i] <= 1); // interdiction constraint

				IloExpr cons(env);
				for (size_t j = 0; j < neighborhood_[i].size(); j++) cons += yVars[neighborhood_[i][j]];

				neighbor_cons.add(cons >= 1);
				cons.end();
			}
			model.add(neighbor_cons);


			// S_k constraints

			// mixing constraint
			for (int i = 0; i < nnode_; i++)
			{
				IloExpr mixing_cons(env);
				int mu_i = neighborhood_[i].size();
				for (size_t j = 0; j < neighborhood_[i].size(); j++)	mixing_cons += yVars[neighborhood_[i][j]];

				for (size_t t = 0; t < rhweight[i].size(); t++)
				{
					model.add(mixing_cons + (rhweight[i][t] - mu_i) * zVars[i][t] <= rhweight[i][t]);
				}
				mixing_cons.end();
			}

			// logitstic constraints
			for (size_t t = 0; t < left_T1.size(); t++)
			{
				IloExpr lower_cons(env), upper_cons(env);
				vector<int> T_1 = left_T1[t], T_2 = right_T2[t];

				for (size_t j = 0; j < T_1.size(); j++) upper_cons -= (xVars[T_1[j]] + yVars[T_1[j]]);
				for (size_t j = 0; j < T_2.size(); j++) upper_cons -= (1 - yVars[T_2[j]] + xVars[T_2[j]]);

				for (int i = 0; i < nnode_; i++)
				{
					if (hweight[i][index_t[i][t]] <= 0)
					{
						lower_cons += 1;
						upper_cons += 1;
					}
					else
					{
						lower_cons += zVars[i][index_rh[i][index_t[i][t]]];
						upper_cons += zVars[i][index_rh[i][index_t[i][t]]];
					}
				}

				//model.add(lower_cons >= nnode_ - 1);
				model.add(upper_cons <= nnode_ - 1);
				lower_cons.end();
				upper_cons.end();
			}



			// cplex solver
			IloCplex cplex(model);

			cplex.setParam(IloCplex::TiLim, 1800);
			//cplex.setParam(IloCplex::Param::Threads, 1);
			string cplex_log = log_path_ + "/cplex_log_" + filename_ + "_" + to_string(K) + ".log";
			ofstream out_file(cplex_log, std::ios_base::out);
			cplex.setOut(out_file);

			cplex.solve();
			

			IloNum ub = floor(cplex.getObjValue());
			IloNumArray opt_x(env), opt_y(env);
			cplex.getValues(opt_x, xVars);
			cplex.getValues(opt_y, yVars);
			//env.out() << "Optimal x: " << opt_x << endl;
			//env.out() << "Optimal y: " << opt_y << endl;
			
			/**********************************************************
			*
			*				Solve the lower-level problem
			*
			***********************************************************/
			IloModel lower_model(env);
			lower_model.add(yVars);
			// interdiction constraint
			lower_model.add(IloScalProd(opt_x, yVars) == 0);
			lower_model.add(neighbor_cons);
			neighbor_cons.end();

			// objective
			lower_model.add(IloMinimize(env, IloScalProd(follower_weight_, yVars), "lower-level objective function"));

			cplex.extract(lower_model);
			cplex.solve();
			clock_t end_t = clock();

			cplex.getValues(opt_y, yVars);
			//env.out() << "Optimal lower-level y: " << opt_y << endl;

			IloNum lb = IloScalProd(leader_weight_, opt_y);
			//env.out() << "Lower Bound: " << lb << endl;

			//double gap = (ub - lb) / float(ub);
			std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

			// record the result
			cplex.out() << K << '\t' << ub << '\t' << lb << '\t' << float(std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count()) / 1000 << endl;
			cout << K << '\t' << ub << '\t' << lb << '\t' << float(std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count()) / 1000 << endl;
			out_result << ub << '\t' << lb << '\t' << float(std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count()) / 1000 << '\t';

			// cplex.out() << K << '\t' << ub << '\t' << lb << '\t' << float(end_t - start_t) / CLOCKS_PER_SEC << endl;
			// cout << K << '\t' << ub << '\t' << lb << '\t' << float(end_t - start_t) / CLOCKS_PER_SEC << endl;
			// out_result << ub << '\t' << lb << '\t' << float(end_t - start_t) / CLOCKS_PER_SEC << '\t';

			// clear
			cplex.end();
			lower_model.end();
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
	bool Solver_local_ext(int K, ofstream &out_result) {
		/**********************************************************
		*
		*				Find T^k
		*
		***********************************************************/
		vector<vector<int> > left_T1, right_T2;
		vector<vector<int> > hweight;

		vector<vector<size_t> > index_h, index_t;
		vector<vector<int> > rhweight(nnode_), index_rh;


		clock_t start_t = clock();
		std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
		// find T for K>=2 
		if (K >= 1)
		{
			local_pair_k(left_T1, right_T2, hweight, K);

			for (size_t i = 0; i < nnode_; i++)
			{
				vector<int> hweight_i = hweight[i];
				vector<size_t> idx_i(hweight[i].size()), idx_ti(hweight[i].size());
				iota(idx_i.begin(), idx_i.end(), 0);
				iota(idx_ti.begin(), idx_ti.end(), 0);

				// sort indexes based on comparing values in v
				sort(idx_i.begin(), idx_i.end(),
					[&hweight_i](size_t i1, size_t i2) {return hweight_i[i1] < hweight_i[i2]; });
				sort(idx_ti.begin(), idx_ti.end(),
					[&idx_i](size_t i1, size_t i2) {return idx_i[i1] < idx_i[i2]; });

				sort(hweight[i].begin(), hweight[i].end());
				index_h.push_back(idx_i);
				index_t.push_back(idx_ti);

				// reduce the weight size
				vector<int> rhweight_i; // remove the duplication items in hweight
				vector<int> idx_ri(hweight[i].size(), -1); // the index that connect hweight to rhweight
				for (size_t t = 0; t < hweight[i].size(); t++) {
					if (hweight[i][t] <= 0) continue;

					if (t == 0) {
						rhweight_i.push_back(hweight[i][t]);
						idx_ri[0] = 0;
						continue;
					}

					if (hweight[i][t] != hweight[i][t - 1])	rhweight_i.push_back(hweight[i][t]);

					idx_ri[t] = rhweight_i.size() - 1;
				}

				index_rh.push_back(idx_ri);
				rhweight[i] = (rhweight_i);
			}
		}


		/**********************************************************
		*
		*				Modeling
		*
		***********************************************************/


		try
		{
			IloModel model(env);

			// create varaibles
			IloNumVarArray xVars(env);
			IloNumVarArray yVars(env);
			NumVarMatrix zVars(env, nnode_);
			NumVarMatrix vVars(env, nnode_);
			
			for (size_t i = 0; i < nnode_; i++)
			{
				xVars.add(IloNumVar(env, 0, 1, IloNumVar::Bool, create_name("x", i).c_str()));
				yVars.add(IloNumVar(env, 0, 1, IloNumVar::Bool, create_name("y", i).c_str()));

				zVars[i] = IloNumVarArray(env);
				vVars[i] = IloNumVarArray(env);
				for (size_t t = 0; t < rhweight[i].size(); t++)
				{
					//zVars[i].add(IloNumVar(env, 0, 1, IloNumVar::Bool, create_name("z", i, t).c_str()));
					zVars[i].add(IloNumVar(env, 0, 1, IloNumVar::Bool, create_name("z", i, t).c_str()));
					vVars[i].add(IloNumVar(env, 0, 1, IloNumVar::Bool, create_name("v", i, t).c_str()));
				}
				model.add(zVars[i]);
				model.add(vVars[i]);
			}

			model.add(xVars);
			model.add(yVars);

			// objective function
			model.add(IloMaximize(env, IloScalProd(leader_weight_, yVars), "Objective"));

			// S_0 constraints
			// leader budget
			model.add(IloSum(xVars) <= budget_);
			IloRangeArray neighbor_cons(env);
			for (size_t i = 0; i < nnode_; i++)
			{
				model.add(xVars[i] + yVars[i] <= 1); // interdiction constraint

				IloExpr cons(env);
				for (size_t j = 0; j < neighborhood_[i].size(); j++) cons += yVars[neighborhood_[i][j]];

				neighbor_cons.add(cons >= 1);
				cons.end();
			}
			model.add(neighbor_cons);


			// S_k constraints

			// mixing constraint
			for (int i = 0; i < nnode_; i++)
			{
				for (size_t t = 0; t < rhweight[i].size(); t++)
				{
					// precedence constraint
					model.add(vVars[i][t] - zVars[i][t] <= 0);
					if (t < rhweight[i].size() - 1)	model.add(vVars[i][t] - vVars[i][t + 1] >= 0);
				
					//if (t < rhweight[i].size() - 1)	model.add(zVars[i][t] - zVars[i][t + 1] >= 0);
				}


				// mixing constraint
				IloExpr mixing_cons(env);
				int mu_i = neighborhood_[i].size();
				for (size_t j = 0; j < neighborhood_[i].size(); j++)	mixing_cons += yVars[neighborhood_[i][j]];

				for (size_t t = 0; t < rhweight[i].size() - 1; t++)
				{
					mixing_cons += (rhweight[i][t] - rhweight[i][t + 1]) * vVars[i][t];
					//mixing_cons += (rhweight[i][t] - rhweight[i][t + 1]) * zVars[i][t];
				}
				mixing_cons += (rhweight[i][rhweight[i].size() - 1] - mu_i) * vVars[i][rhweight[i].size() - 1];


				//mixing_cons += (rhweight[i][rhweight[i].size() - 1] - mu_i) * zVars[i][rhweight[i].size() - 1];

				model.add(mixing_cons <= rhweight[i][0]);
				mixing_cons.end();
			}

			// logitstic constraints
			for (size_t t = 0; t < left_T1.size(); t++)
			{
				IloExpr lower_cons(env), upper_cons(env);
				vector<int> T_1 = left_T1[t], T_2 = right_T2[t];

				for (size_t j = 0; j < T_1.size(); j++) upper_cons -= (xVars[T_1[j]] + yVars[T_1[j]]);
				for (size_t j = 0; j < T_2.size(); j++) upper_cons -= (1 - yVars[T_2[j]] + xVars[T_2[j]]);

				for (int i = 0; i < nnode_; i++)
				{
					if (hweight[i][index_t[i][t]] <= 0)
					{
						lower_cons += 1;
						upper_cons += 1;
					}
					else
					{
						lower_cons += zVars[i][index_rh[i][index_t[i][t]]];
						upper_cons += zVars[i][index_rh[i][index_t[i][t]]];
					}
				}

				//model.add(lower_cons >= nnode_ - 1);
				model.add(upper_cons <= nnode_ - 1);
				lower_cons.end();
				upper_cons.end();
			}



			// cplex solver
			IloCplex cplex(model);
			cplex.setParam(IloCplex::TiLim, 1800);
			//cplex.setParam(IloCplex::Param::Threads, 1);
			string cplex_log = log_path_ + "/cplex_log_" + filename_ + "_" + to_string(K) + "_ext.log";
			ofstream out_file(cplex_log, std::ios_base::out);
			cplex.setOut(out_file);
			
			cplex.solve();
			

			IloNum ub = floor(cplex.getObjValue());
			//env.out() << "Upper Bound: " << ub << endl;
			IloNumArray opt_x(env), opt_y(env);
			cplex.getValues(opt_x, xVars);
			cplex.getValues(opt_y, yVars);
			/*env.out() << "Optimal x: " << opt_x << endl;
			env.out() << "Optimal y: " << opt_y << endl;*/
			//env.out() << "Optimal w: " << opt_w << endl;

			/**********************************************************
			*
			*				Solve the lower-level problem
			*
			***********************************************************/
			IloModel lower_model(env);
			lower_model.add(yVars);
			// interdiction constraint
			lower_model.add(IloScalProd(opt_x, yVars) == 0);
			lower_model.add(neighbor_cons);
			neighbor_cons.end();

			// objective
			lower_model.add(IloMinimize(env, IloScalProd(follower_weight_, yVars), "lower-level objective function"));

			cplex.extract(lower_model);
			cplex.solve();

			clock_t end_t = clock();

			cplex.getValues(opt_y, yVars);
			//env.out() << "Optimal lower-level y: " << opt_y << endl;

			IloNum lb = IloScalProd(leader_weight_, opt_y);
			//env.out() << "Lower Bound: " << lb << endl;

			//double gap = (ub - lb) / float(ub);
			std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

			// record the result
			cplex.out() << K << '\t' << ub << '\t' << lb << '\t' << float(std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count()) / 1000 << endl;
			cout << K << '\t' << ub << '\t' << lb << '\t' << float(std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count()) / 1000 << endl;
			out_result << float(std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count()) / 1000 << '\t';

			// // record the result
			// cplex.out() << K << '\t' << ub << '\t' << lb << '\t' << float(end_t - start_t) / CLOCKS_PER_SEC << endl;
			// cout << K << '\t' << ub << '\t' << lb << '\t' << float(end_t - start_t) / CLOCKS_PER_SEC << endl;
			// out_result << float(end_t - start_t) / CLOCKS_PER_SEC << '\t';
			// clear
			cplex.end();
			lower_model.end();
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
	void local_pair_k(vector<vector<int> > &left_T1, vector<vector<int> > &right_T2, vector<vector<int> > &hweight, int K) {
		vector<int> N;
		hweight.resize(nnode_);

		// K = 1
		for (size_t i = 0; i < nnode_; i++)
		{
			N.push_back(i);

			left_T1.push_back(vector<int>());
			right_T2.push_back(vector<int>(1, i));

			for (size_t j = 0; j < nnode_; j++)
			{
				if (find(neighborhood_[j].begin(), neighborhood_[j].end(), i) != neighborhood_[j].end())
				{
					hweight[j].push_back(1);
				}
				else
				{
					hweight[j].push_back(0);
				}
			}
		}

		// K = 2, ..., K
		for (int k = 2; k <= K; k++)
		{
			// t is the size of one set
			for (size_t t = 1; t <= floor(float(k) / 2.0); t++)
			{
				vector<vector<int> > lcomb = subset_k(N, t);

				for (size_t liter = 0; liter < lcomb.size(); liter++)
				{
					vector<int> T_1 = lcomb[liter];

					int wsum_T1(0), min_T1 = std::numeric_limits<int>::max();
					for (size_t i = 0; i < T_1.size(); i++)
					{
						wsum_T1 += follower_weight_[T_1[i]];
						if (min_T1 > follower_weight_[T_1[i]])	min_T1 = follower_weight_[T_1[i]];
					}

					vector<int> N_T; // the complement of T1
					set_difference(N.begin(), N.end(), T_1.begin(), T_1.end(),
						std::inserter(N_T, N_T.begin()));


					//printVec(T_1);
					//cout << wsum_T1 << '\t' << min_T1 << endl;


					vector<vector<int> > rcomb = subset_k(N_T, k - t);
					for (size_t riter = 0; riter < rcomb.size(); riter++)
					{
						vector<int> T_2 = rcomb[riter];

						int wsum_T2(0), min_T2 = std::numeric_limits<int>::max();
						for (size_t i = 0; i < T_2.size(); i++)
						{
							wsum_T2 += follower_weight_[T_2[i]];
							if (min_T2 > follower_weight_[T_2[i]])	min_T2 = follower_weight_[T_2[i]];
						}


						//printVec(T_2);
						//cout << wsum_T2 << '\t' << min_T2 << endl;

						// remove if N_i \subseteq T1 or T2
						bool flag = 0;
						vector<int> hweight_t(nnode_);
						for (size_t i = 0; i < nnode_; i++)
						{
							vector<int> inter_v1, inter_v2;
							set_intersection(
								T_1.begin(), T_1.end(), neighborhood_[i].begin(), neighborhood_[i].end(),
								inserter(inter_v1, inter_v1.begin()));
							set_intersection(
								T_2.begin(), T_2.end(), neighborhood_[i].begin(), neighborhood_[i].end(),
								inserter(inter_v2, inter_v2.begin()));

							if ((inter_v1.size() == neighborhood_[i].size()) || (inter_v2.size() == neighborhood_[i].size())) {
								flag = 1;
								break;
							}

							hweight_t[i] = inter_v2.size() - inter_v1.size();
						}
						if (flag) continue;


						// determine left and right
						if ((wsum_T1 < wsum_T2) && (wsum_T2 - min_T2 <= wsum_T1))
						{
							left_T1.push_back(T_1);
							right_T2.push_back(T_2);

							for (size_t i = 0; i < nnode_; i++)	hweight[i].push_back(hweight_t[i]);
						}
						else if ((float(t) != float(k) / 2.0) && (wsum_T1 > wsum_T2) && (wsum_T1 - min_T1 <= wsum_T2))
						{
							left_T1.push_back(T_2);
							right_T2.push_back(T_1);

							for (size_t i = 0; i < nnode_; i++)	hweight[i].push_back(-hweight_t[i]);
						}

					}

				}

			}

		}
	}

	string filename_;
	string log_path_;
	int budget_;
	int nnode_;
	int nedge_;

	vector<vector<int> > vedge_;
	vector<vector<int> > neighborhood_;

	IloEnv env;
	IloIntArray leader_weight_;
	IloIntArray follower_weight_;
};


int main() {


	/********************************************************************************
	*
	*                               Bilevel Vertex Cover
	*
	*********************************************************************************/

	/*
	* 	test for one instance
	*/
	// string file_path = "./Dataset/BVC/";
	// string file_name = "BVC_n15d8b8j7";
	// string log_path = "Log/Log_BVC_ASymm";
	// BVC test_prob(file_name, file_path, log_path, 0);

	// ofstream out_file("results.txt");
	
	// test_prob.Solver_local(3, out_file);
	// test_prob.Solver_local_ext(3, out_file);
	//test_prob.Solver_CP(out_file);

	string file_path = "./Dataset/BVC/";
	string log_path;


	// Experiments for n=15
	for (int flag = 0; flag <= 1; ++flag)
	{
		string outfile_name;
		if(flag){
			outfile_name = "BVC_n20_Symm_results.txt";
			log_path = "Log/Log_BVC_Symm";
		}else{
			outfile_name = "BVC_n20_Asymm_results.txt";
			log_path = "Log/Log_BVC_ASymm";
		}

		ofstream out_file(outfile_name);
		out_file << "\t | \t \t BVC_0 \t | \t \t BVC_1 \t \t | \t \t BVC_2 \t \t | \t \t BVC_3 \t \t |" << endl;
		out_file << "File | \t ObjU \t ObjL \t Time | \t ObjU \t ObjL \t Time \t ExtTime | \t ObjU \t ObjL \t Time \t ExtTime | \t ObjU \t ObjL \t Time \t ExtTime |" << endl;

		for (int i = 0; i < 10; ++i)
		{
			int n = 20;
			int deg = ceil(n / float(2));

			string file_name = "BVC_n" + to_string(n) + "d" + to_string(deg) +  "b" + to_string(deg) + "j" + to_string(i);

			cout << file_name << endl;
			BVC test_prob(file_name, file_path, log_path, flag);
			out_file << file_name << "\t";

			for (int k = 0; k < 8; k++)
			{
				test_prob.Solver_local(k, out_file);
			}
	   
	    	out_file << '\n';
		}

		out_file.close();

	}



	// // Experiments for all test instances
	// for (int flag = 0; flag <= 1; ++flag)
	// {
	// 	string outfile_name;
	// 	if(flag){
	// 		outfile_name = "BVC_Symm_results.txt";
	// 		log_path = "Log/Log_BVC_Symm";
	// 	}else{
	// 		outfile_name = "BVC_Asymm_results.txt";
	// 		log_path = "Log/Log_BVC_ASymm";
	// 	}

	// 	ofstream out_file(outfile_name);
	// 	out_file << "\t | \t \t BVC_0 \t | \t \t BVC_1 \t \t | \t \t BVC_2 \t \t | \t \t BVC_3 \t \t |" << endl;
	// 	out_file << "File | \t ObjU \t ObjL \t Time | \t ObjU \t ObjL \t Time \t ExtTime | \t ObjU \t ObjL \t Time \t ExtTime | \t ObjU \t ObjL \t Time \t ExtTime |" << endl;

	// 	for (int n = 10; n <= 50; n += 5)
	// 	{
	// 		int deg = ceil(n / float(2));
	// 		for (int i = 0; i < 10; i++)
	// 		{
	// 			string file_name = "BVC_n" + to_string(n) + "d" + to_string(deg) +  "b" + to_string(deg) + "j" + to_string(i);

	// 			cout << file_name << endl;
	// 			BVC test_prob(file_name, file_path, log_path, flag);

	// 			out_file << file_name << "\t";
	// 			test_prob.Solver_local(0, out_file);
	// 			test_prob.Solver_local(1, out_file);
	// 			test_prob.Solver_local_ext(1, out_file);
	// 			test_prob.Solver_local(2, out_file);
	// 			test_prob.Solver_local_ext(2, out_file);
	// 			test_prob.Solver_local(3, out_file);
	// 			test_prob.Solver_local_ext(3, out_file);
	// 			out_file << '\n';

	// 		}
	// 	}

	// 	out_file.close();
	// }


	// // Experiments for n=20, b=10, deg in [10, 15]
	// for (int flag = 0; flag <= 1; ++flag)
	// {
	// 	string outfile_name;
	// 	if(flag){
	// 		outfile_name = "BVC_n20_Symm_results.txt";
	// 		log_path = "Log/Log_BVC_Symm";
	// 	}else{
	// 		outfile_name = "BVC_n20_Asymm_results.txt";
	// 		log_path = "Log/Log_BVC_ASymm";
	// 	}

	// 	ofstream out_file(outfile_name);
	// 	out_file << "\t | \t \t BVC_0 \t | \t \t BVC_1 \t \t | \t \t BVC_2 \t \t | \t \t BVC_3 \t \t |" << endl;
	// 	out_file << "File | \t ObjU \t ObjL \t Time | \t ObjU \t ObjL \t Time \t ExtTime | \t ObjU \t ObjL \t Time \t ExtTime | \t ObjU \t ObjL \t Time \t ExtTime |" << endl;

	// 	int n = 20, b = 10;
	// 	for (int deg = 10; deg <= 15; deg++)
	// 	{
	// 		for (int i = 0; i < 10; i++)
	// 		{
	// 			string file_name = "BVC_n" + to_string(n) + "d" + to_string(deg) +  "b" + to_string(b) + "j" + to_string(i);

	// 			cout << file_name << endl;
	// 			BVC test_prob(file_name, file_path, log_path, flag);

	// 			out_file << file_name << "\t";
	// 			test_prob.Solver_local(0, out_file);
	// 			test_prob.Solver_local(1, out_file);
	// 			test_prob.Solver_local_ext(1, out_file);
	// 			test_prob.Solver_local(2, out_file);
	// 			test_prob.Solver_local_ext(2, out_file);
	// 			test_prob.Solver_local(3, out_file);
	// 			test_prob.Solver_local_ext(3, out_file);
	// 			out_file << '\n';

	// 		}
	// 	}

	// 	out_file.close();
	// }


	return 0;
}