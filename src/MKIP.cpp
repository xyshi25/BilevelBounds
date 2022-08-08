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

typedef IloArray<IloIntArray> IntMatrix;

typedef IloArray<IloNumVarArray> IloVarMatrix;

// // print vector
// template<typename T>
// void printVec(const vector<T> &vec) {

// 	cout << "[ ";
// 	for (size_t i = 0; i < vec.size(); i++)
// 	{
// 		cout << vec[i] << " ";
// 	}

// 	cout << "]" << endl;
// }

// // Print matrix
// template<typename T>
// void printMat(const vector<vector<T> > &vec) {

// 	cout << "[ " << endl;
// 	for (size_t i = 0; i < vec.size(); i++)
// 	{
// 		cout << i << ": ";
// 		for (size_t j = 0; j < vec[i].size(); j++)
// 		{
// 			cout << vec[i][j] << " ";
// 		}
// 		cout << endl;
// 	}

// 	cout << "]" << endl;
// }

class BMKP
{
public:
	BMKP()
		: n_(0), m1_(0), m2_(0), ubudget_(env), lbudget_(env), mu_(env), cost_(env), uweight_(env), lweight_(env) {}
	
	BMKP(string filename, string file_path, string log_path) {
		log_path_ = log_path;
		ubudget_ = IloIntArray(env);
		lbudget_ = IloIntArray(env);
		mu_ = IloIntArray(env);

		cost_ = IloIntArray(env);
		uweight_ = IntMatrix(env);
		lweight_ = IntMatrix(env);
		if (!ReadFile(filename, file_path)) {
			return;
		}
	}
	~BMKP() {
		Clear();
	}

	void Clear() {
		n_ = 0;
		m1_ = 0;
		m2_ = 0;

		ubudget_.clear(); 
		lbudget_.clear();
		cost_.clear();

		for (size_t i = 0; i < m1_; i++) uweight_.clear();
		uweight_.clear();

		for (size_t i = 0; i < m2_; i++) lweight_.clear();
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
					split >> n_ >> m1_ >> m2_;
					continue;
				}

				// ubudget_

				// lbudget_

				vector<string> split_line((istream_iterator<string>(split)), istream_iterator<string>());

				if (split_line.size() == 0) continue;

				if (split_line[0] == "ubudget")
				{
					for (size_t i = 1; i < split_line.size(); i++) ubudget_.add(IloInt(stoi(split_line[i])));
					assert(ubudget_.getSize() == m1_);
				}
				else if (split_line[0] == "lbudget")
				{
					for (size_t i = 1; i < split_line.size(); i++) lbudget_.add(IloInt(stoi(split_line[i])));
					assert(lbudget_.getSize() == m1_);
				}
				else if (split_line[0] == "uweight")
				{
					for (int i = 0; i < m1_; ++i)
					{
						getline(in_file, line);
						istringstream split(line);
						line_num++;
						
						vector<string> split_line((istream_iterator<string>(split)), istream_iterator<string>());

						IloIntArray weight(env);
						for (size_t i = 0; i < split_line.size(); i++) weight.add(IloInt(stoi(split_line[i])));
						assert(lbudget_.getSize() == n_);
						uweight_.add(weight);
					}
				}
				else if (split_line[0] == "lweight")
				{
					for (int i = 0; i < m2_; ++i)
					{
						getline(in_file, line);
						istringstream split(line);
						line_num++;
						
						vector<string> split_line((istream_iterator<string>(split)), istream_iterator<string>());

						IloIntArray weight(env);
						for (size_t i = 0; i < split_line.size(); i++) weight.add(IloInt(stoi(split_line[i])));
						assert(lbudget_.getSize() == n_);
						lweight_.add(weight);
					}
				}
				else if (split_line[0] == "cost")
				{
					for (size_t i = 1; i < split_line.size(); i++) cost_.add(stoi(split_line[i]));
					assert(cost_.getSize() == n_);
				}

			}

			//env.out() << n_ << '\t' << m1_ << '\t' << m2_ << '\n' << cost_ << '\n' << ubudget_ << '\n' << lbudget_ << '\n' << uweight_ << '\n' << lweight_ << '\n';

			// compute mu
			for (int i = 0; i < lweight_.getSize(); ++i)
			{
				mu_.add(0);
			}
			
			return 1;
		}
		else
		{
			cout << "Can't open file " << filename << endl;
			return 0;
		}
	}
	
	bool WriteToMps(string filename) {
		string mps_file = "Dataset/MKIPMibs/" + filename + ".mps";
		string aux_file = "Dataset/MKIPMibs/" + filename + ".txt";
		ofstream mps_out(mps_file);


		mps_out << "NAME\t" << filename << "\n";
		//ROWS
		mps_out << "ROWS\n N \t OBJROW\n";
		mps_out.fill('0');
		for (size_t i = 0; i < m2_; i++)
		{
			mps_out << " L\tR" << setw(4) << i << '\n';
		}

		// COLUMNS
		mps_out << "COLUMNS\n";

		for (size_t i = 0; i < n_; i++)
		{
			mps_out << " C" << setw(4) << i << '\t';
			mps_out << "OBJROW\t" << -cost_[i] << '\t';
			for (size_t j = 0; j < m2_; j++)
			{
				mps_out << "R" << setw(4) << j << '\t' << lweight_[j][i] << '\t';
			}
			mps_out << endl;
		}

		// RHS
		mps_out << "RHS\n";
		for (size_t i = 0; i < m2_; i++)
		{
			mps_out << " RHS\t R" << setw(4) << i << '\t' << lbudget_[i] << endl;
		}

		// BOUNDS
		mps_out << "BOUNDS\n";

		for (size_t i = 0; i < n_; i++)
		{
			mps_out << " BV BOUND\t" << " C" << setw(4) << i << "\t 1" << "\t\n";
		}

		mps_out << "ENDATA\n";
		mps_out.close();
		///////////////////////////////// END mps_out

		ofstream aux_out(aux_file);
		aux_out << "N " << n_ << "\nM " << n_ + m2_ << endl;
		for (size_t i = 0; i < n_; i++)	aux_out << "LC " << n_ + i << endl;

		for (size_t i = 0; i < n_+m2_; i++) aux_out << "LR " << 1 + i << endl;

		for (size_t i = 0; i < n_; i++) aux_out << "LO " << -cost_[i] << endl;
		aux_out << "OS " << 1 << endl;

		for (size_t i = 0; i < n_; i++) aux_out << "IC " << uweight_[0][i] << endl;

		aux_out << "IB " << ubudget_[0] << endl;
		aux_out.close();


		return 1;
	}
	
	
	
	bool Solver_local_ext(int K, ofstream &out_result) {
		clock_t start_t = clock();
		std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

		vector<vector<int> > left_T, right_T;
		vector<vector<int> > hweight;
		vector<vector<int> > rhweight(m2_);
		vector<vector<size_t> > idx_h, idx_t, idx_rh;
		// find T for K>=2 
		if (K >= 1)
		{
			local_pair_k(left_T, right_T, hweight, K);

			for (int i = 0; i < m2_; ++i)
			{
				vector<int> value = hweight[i];

				vector<size_t> idx_v(value.size()), idx_tv(value.size());
				iota(idx_v.begin(), idx_v.end(), 0);

				// sort indexes based on comparing values in v
				sort(idx_v.begin(), idx_v.end(),
					[&value](size_t i1, size_t i2) {return value[i1] < value[i2]; });

				for (int j = 0; j < value.size(); ++j)	idx_tv[idx_v[j]] = j;
				idx_h.push_back(idx_v);
				idx_t.push_back(idx_tv);

				std::vector<int> rvalue;
				std::vector<size_t> idx_rv(value.size());
				for (int j = 0; j < value.size(); ++j)
				{
					if (value[idx_v[j]] <= 0)
					{
						idx_rv[j + 1] = 0;
						continue;
					}

					if (j == 0 || (value[idx_v[j]] != value[idx_v[j - 1]]))
					{
						rvalue.push_back(value[idx_v[j]]);
					}

					idx_rv[j] = rvalue.size() - 1;
				}

				rhweight[i] = (rvalue);
				idx_rh.push_back(idx_rv);
			}
		}

		try
		{
			IloModel model(env);

			//  create varaibles
			IloNumVarArray xVars(env);
			IloNumVarArray yVars(env);
			IloVarMatrix wVars(env, m2_);
			for (size_t i = 0; i < n_; i++)
			{
				xVars.add(IloNumVar(env, 0, 1, IloNumVar::Bool, create_name("x", i).c_str()));
				yVars.add(IloNumVar(env, 0, 1, IloNumVar::Bool, create_name("y", i).c_str()));
			}
			model.add(xVars);
			model.add(yVars);


			for (int i = 0; i < m2_; ++i)
			{
				wVars[i] = IloNumVarArray(env);
				for (int j = 0; j < rhweight[i].size(); ++j)
				{
					wVars[i].add(IloNumVar(env, 0, 1, IloNumVar::Bool, create_name("w", i, j).c_str()));
				}
				model.add(wVars[i]);

			}
			

			// objective function
			model.add(IloMinimize(env, IloScalProd(cost_, yVars), "objective function"));

			// S_0 constraints
			for (size_t i = 0; i < n_; i++)
			{
				model.add(xVars[i] + yVars[i] <= 1);
			}

			for (int i = 0; i < m1_; ++i)
			{
				model.add(IloScalProd(uweight_[i], xVars) <= ubudget_[i]);
			}

			for (int i = 0; i < m2_; ++i)
			{
				model.add(IloScalProd(lweight_[i], yVars) <= lbudget_[i]);
			}

			// S_k constraints
			
			for (size_t i = 0; i < m2_; i++)
			{
				IloExpr mixing_con(env);
				mixing_con = IloScalProd(lweight_[i], yVars);
				for (size_t j = 0; j < rhweight[i].size(); j++)
				{
					if (j == 0)
					{
						mixing_con -= rhweight[i][j] * wVars[i][j];
					}
					else
					{
						mixing_con += (rhweight[i][j - 1] - rhweight[i][j]) * wVars[i][j];
					}

					if (j < rhweight[i].size() - 1)	model.add(wVars[i][j] - wVars[i][j + 1] >= 0);

				}
				model.add(mixing_con >= 0);
				mixing_con.end();
			}


			for (size_t t = 0; t < left_T.size(); t++)
			{
				vector<int> T_1 = left_T[t], T_2 = right_T[t];
				IloExpr upper_cons(env);

				for (int i = 0; i < T_1.size(); ++i)  upper_cons += xVars[T_1[i]] + yVars[T_1[i]];
				for (int i = 0; i < T_2.size(); ++i) upper_cons += xVars[T_2[i]] + 1 - yVars[T_2[i]];

				for (int i = 0; i < m2_; ++i)
				{
					if (hweight[i][t] <= 0) continue;

					upper_cons += wVars[i][idx_rh[i][idx_t[i][t]]];
				}
				model.add(upper_cons >= 1);
				upper_cons.end();
			}

			// cplex solver
			IloCplex cplex(model);
			cplex.setParam(IloCplex::TiLim, 600);
			//cplex.setParam(IloCplex::Param::Threads, 1);
			string cplex_log = log_path_ + "/cplex_log_" + filename_ + "_" + to_string(K) + "_ext.log";
			ofstream out_file(cplex_log, std::ios_base::out);
			cplex.setOut(out_file);
			
			cplex.solve();
			clock_t end_t = clock();;

			IloNum lb = cplex.getBestObjValue();
			IloNumArray opt_x(env), opt_y(env);
			cplex.getValues(opt_x, xVars);
			//cplex.getValues(opt_y, yVars);

			// solve the lower level problem
			IloModel lower_model(env);
			lower_model.add(yVars);

			for (int i = 0; i < m2_; ++i)
			{
				lower_model.add(IloScalProd(lweight_[i], yVars) <= lbudget_[i]);
			}

			// interdiction constraint
			lower_model.add(IloScalProd(opt_x, yVars) == 0);

			// objective for lower level
			lower_model.add(IloMaximize(env, IloScalProd(cost_, yVars), "lower-level objective function"));

			cplex.extract(lower_model);
			cplex.solve();

			IloNum ub = cplex.getObjValue();
			
			double gap = (ub - lb) / float(ub);

			std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
			cout << K << '\t' << lb << '\t' << ub << '\t' << float(std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count()) / 1000 << endl;
			out_result  << left_T.size() << '\t' << lb << '\t' << ub << '\t' << float(std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count()) / 1000 << '\t';

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
		clock_t start_t = clock();
		std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

		vector<vector<int> > left_T, right_T;
		vector<vector<int> > hweight;
		vector<vector<int> > rhweight(m2_);
		vector<vector<size_t> > idx_h, idx_t, idx_rh;
		// find T for K>=2 
		if (K >= 1)
		{
			local_pair_k(left_T, right_T, hweight, K);

			for (int i = 0; i < m2_; ++i)
			{
				vector<int> value = hweight[i];

				vector<size_t> idx_v(value.size()), idx_tv(value.size());
				iota(idx_v.begin(), idx_v.end(), 0);

				// sort indexes based on comparing values in v
				sort(idx_v.begin(), idx_v.end(),
					[&value](size_t i1, size_t i2) {return value[i1] < value[i2]; });

				for (int j = 0; j < value.size(); ++j)	idx_tv[idx_v[j]] = j;
				idx_h.push_back(idx_v);
				idx_t.push_back(idx_tv);

				std::vector<int> rvalue;
				std::vector<size_t> idx_rv(value.size());
				for (int j = 0; j < value.size(); ++j)
				{
					if (value[idx_v[j]] <= 0)
					{
						idx_rv[j + 1] = 0;
						continue;
					}

					if (j == 0 || (value[idx_v[j]] != value[idx_v[j - 1]]))
					{
						rvalue.push_back(value[idx_v[j]]);
					}

					idx_rv[j] = rvalue.size() - 1;
				}

				rhweight[i] = (rvalue);
				idx_rh.push_back(idx_rv);
			}
		}


		try
		{
			IloModel model(env);

			//  create varaibles
			IloNumVarArray xVars(env);
			IloNumVarArray yVars(env);
			IloVarMatrix wVars(env, m2_);
			for (size_t i = 0; i < n_; i++)
			{
				xVars.add(IloNumVar(env, 0, 1, IloNumVar::Bool, create_name("x", i).c_str()));
				yVars.add(IloNumVar(env, 0, 1, IloNumVar::Bool, create_name("y", i).c_str()));
			}
			model.add(xVars);
			model.add(yVars);


			for (int i = 0; i < m2_; ++i)
			{
				wVars[i] = IloNumVarArray(env);
				for (int j = 0; j < rhweight[i].size(); ++j)
				{
					wVars[i].add(IloNumVar(env, 0, 1, IloNumVar::Bool, create_name("w", i, j).c_str()));
				}
				model.add(wVars[i]);
				
			}


			// objective function
			model.add(IloMinimize(env, IloScalProd(cost_, yVars), "objective function"));

			// S_0 constraints
			for (size_t i = 0; i < n_; i++)
			{
				model.add(xVars[i] + yVars[i] <= 1);
			}

			for (int i = 0; i < m1_; ++i)
			{
				model.add(IloScalProd(uweight_[i], xVars) <= ubudget_[i]);
			}

			for (int i = 0; i < m2_; ++i)
			{
				model.add(IloScalProd(lweight_[i], yVars) <= lbudget_[i]);
			}


			// S_k constraints
			for (int i = 0; i < m2_; ++i)
			{
				for (int j = 0; j < rhweight[i].size(); j++)
				{
					model.add(IloScalProd(lweight_[i], yVars) - wVars[i][j]*rhweight[i][j] >= 0);
				}
			}


			for (size_t t = 0; t < left_T.size(); t++)
			{
				vector<int> T_1 = left_T[t], T_2 = right_T[t];
				IloExpr upper_cons(env);

				for (int i = 0; i < T_1.size(); ++i)  upper_cons += xVars[T_1[i]] + yVars[T_1[i]];
				for (int i = 0; i < T_2.size(); ++i) upper_cons += xVars[T_2[i]] + 1 - yVars[T_2[i]];

				for (int i = 0; i < m2_; ++i)
				{
					if (hweight[i][t] <= 0) continue;
					
					upper_cons += wVars[i][idx_rh[i][idx_t[i][t]]];
				}
				model.add(upper_cons >= 1);
				upper_cons.end();
			}



			// cplex solver
			IloCplex cplex(model);
			
			cplex.setParam(IloCplex::TiLim, 600);
			//cplex.setParam(IloCplex::Param::Threads, 1);
			string cplex_log = log_path_ + "/cplex_log_" + filename_ + "_" + to_string(K) + ".log";
			ofstream out_file(cplex_log, std::ios_base::out);
			cplex.setOut(out_file);
	
			cplex.solve();
			clock_t end_t = clock();

			IloNum lb = cplex.getBestObjValue();
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

			for (int i = 0; i < m2_; ++i)
			{
				lower_model.add(IloScalProd(lweight_[i], yVars) <= lbudget_[i]);
			}

			// interdiction constraint
			lower_model.add(IloScalProd(opt_x, yVars) == 0);

			// objective for lower level
			lower_model.add(IloMaximize(env, IloScalProd(cost_, yVars), "lower-level objective function"));

			cplex.extract(lower_model);
			cplex.solve();

			IloNum ub = cplex.getObjValue();

			std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

			out_result << left_T.size()  << '\t' << lb << '\t' << ub << '\t' << float(std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count()) / 1000 << '\t';
			cout << left_T.size() << '\t' <<  lb << '\t' << ub << '\t' << float(std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count()) / 1000 << endl;

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
	void local_pair_k(vector<vector<int> > &left_T1, vector<vector<int> > &right_T1, vector<vector<int> >&hweight, int K) {
		vector<int> N;

		hweight.resize(m2_);
		
		// K = 1
		for (size_t i = 0; i < n_; i++)
		{
			N.push_back(i);

			for (int j = 0; j < m2_; ++j) hweight[j].push_back(lbudget_[j] + 1 - lweight_[j][i]);
			left_T1.push_back(vector<int>(1, i));
			right_T1.push_back(vector<int>());
			
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
					vector<int> wsum_T1(m2_, 0);
					int sum_T1(0), min_T1 = std::numeric_limits<int>::max();
					for (size_t l = 0; l < t; l++)
					{
						
						sum_T1 += cost_[T_1[l]];
						if (min_T1 > cost_[T_1[l]]) min_T1 = cost_[T_1[l]];
						for (int j = 0; j < m2_; ++j)
						{
							wsum_T1[j] += lweight_[j][T_1[l]];
						}
					}


					set_difference(N.begin(), N.end(), T_1.begin(), T_1.end(),
						std::inserter(N_T, N_T.begin()));

					vector<vector<int> > comb_2 = subset_k(N_T, k - t);
					for (size_t j = 0; j < comb_2.size(); j++)
					{
						vector<int> T_2 = comb_2[j];
						vector<int> wsum_T2(m2_, 0);
						int sum_T2(0), min_T2 = std::numeric_limits<int>::max();
						for (size_t l = 0; l < T_2.size(); l++)
						{
							sum_T2 += cost_[T_2[l]];
							if (min_T2 > cost_[T_2[l]]) min_T2 = cost_[T_2[l]];

							for (int p = 0; p < m2_; ++p)	wsum_T2[p] += lweight_[p][T_2[l]];

						}

						// T_1 > T_2
						bool is_valid = false;
						for (int p = 0; p < m2_; ++p)
						{
							if (lbudget_[p] + 1 - abs(wsum_T2[p] - wsum_T1[p]) > 0)
							{
								is_valid = true;
								break;
							}
						}

						if (!is_valid) continue;


						if ((sum_T1 > sum_T2) && (sum_T1 - min_T1 <= sum_T2))
						{
							left_T1.push_back(T_1);
							right_T1.push_back(T_2);
							
							for (int p = 0; p < m2_; ++p) {
								hweight[p].push_back(lbudget_[p] + 1 + wsum_T2[p] - wsum_T1[p]);
							}

						}
						if ((float(t) != float(k) / 2.0) && (sum_T1 < sum_T2) && (sum_T2 - min_T2 <= sum_T1))
						{
							left_T1.push_back(T_2);
							right_T1.push_back(T_1);
							for (int p = 0; p < m2_; ++p) hweight[p].push_back(lbudget_[p]+ 1 + wsum_T1[p] - wsum_T2[p]);

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
	int m1_;
	int m2_;
	IloIntArray ubudget_;
	IloIntArray lbudget_;
	IloIntArray mu_;


	IloIntArray cost_;
	IntMatrix uweight_;
	IntMatrix lweight_;
	IloEnv env;

};


int main() {


	/*********************************************************************************************
	*
	*									DNEG
	*
	**********************************************************************************************/

	/*string file_name="MKIP_n10k4j3_m1_1_m2_1";
	string file_path = "Dataset/MKIP/";
	string log_path = "Log";

	ofstream out_file("results.txt");

	BMKP test_prob(file_name, file_path, log_path);
	test_prob.Solver_local(0, out_file);
	test_prob.Solver_local(1, out_file);
	test_prob.Solver_local_ext(1, out_file);
	test_prob.Solver_local(2, out_file);
	test_prob.Solver_local_ext(2, out_file);
	test_prob.Solver_local(3, out_file);
	test_prob.Solver_local_ext(3, out_file);*/

	/*test_prob.Solver_local(8, out_file);
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

	/*ofstream out_file;
	string file_path = "Dataset/MKIP/";
	
	// Experiments for n=15, r=4
	out_file.open("results.txt");
	out_file << "\t | CCLW | \t \t \t KIP_0 \t \t | \t \t \t KIP_1 \t \t \t |  \t \t \t KIP_2 \t \t \t | \t \t \t KIP_3 \t \t \t |" << endl;
	out_file << "File \t OptObj |\t Time |\t |T| \t #mixing_set \t ObjL \t ObjU \t Time | \t |T| \t #mixing_set \t ObjL \t ObjU \t Time \t ExtTime | \t |T| \t #mixing_set \t ObjL \t ObjU \t Time \t ExtTime | \t |T| \t #mixing_set \t ObjL \t ObjU \t Time \t ExtTime |" << endl;

	for (int i = 0; i < 10; i++)
	{
		string file_name = "DNeg_n15k4j" + to_string(i);
		cout << file_name << endl;
		BMKP test_prob(file_name, file_path, log_path);
		out_file << file_name << "\t";
		
		for (int k = 0; k < 11; k++)
		{
			test_prob.Solver_local(k, out_file);
		}
   
    	out_file << '\n';

	}

	out_file.close();*/



	// // Experiments for n=20, r=4
	// out_file.open("DNeg_n20r4_results.txt");
	// out_file << "\t | CCLW | \t \t \t KIP_0 \t \t | \t \t \t KIP_1 \t \t \t |  \t \t \t KIP_2 \t \t \t | \t \t \t KIP_3 \t \t \t |" << endl;
	// out_file << "File \t OptObj |\t Time |\t |T| \t #mixing_set \t ObjL \t ObjU \t Time | \t |T| \t #mixing_set \t ObjL \t ObjU \t Time \t ExtTime | \t |T| \t #mixing_set \t ObjL \t ObjU \t Time \t ExtTime | \t |T| \t #mixing_set \t ObjL \t ObjU \t Time \t ExtTime |" << endl;

	// for (int i = 0; i < 10; i++)
	// {
	// 	string file_name = "DNeg_n20k4j" + to_string(i);
	// 	cout << file_name << endl;
	// 	BKP test_prob(file_name, file_path, log_path);
	// 	out_file << file_name << "\t";
	// 	test_prob.Solver_CCLW(out_file);

	// 	for (int k = 0; k < 8; k++)
	// 	{
	// 		test_prob.Solver_local(k, out_file);
	// 	}
   
 //    	out_file << '\n';

	// }

	// out_file.close();


	// ifstream in_file("IMKP_small.ins");
	// string line;
	// string file_name;
	// ofstream out_file;
	// string file_path = "Dataset/IMKP/";
	// string log_path = "Log/Log_IMKP";

	// out_file.open("IMKP_results-0.txt");
	// while (getline(in_file, line))
	// {
		
	// 	line.pop_back();
	// 	cout << line << '\t' <<  line + ".dat-0" << endl;

	// 	file_name = line + ".dat-0";
	// 	cout << file_name << endl;
	// 	BMKP test_prob(file_name, file_path, log_path);

	// 	out_file << file_name << "\t";
	// 	test_prob.Solver_local(0, out_file);
	// 	test_prob.Solver_local(1, out_file);
	// 	test_prob.Solver_local_ext(1, out_file);
	// 	test_prob.Solver_local(2, out_file);
	// 	test_prob.Solver_local_ext(2, out_file);
	// 	test_prob.Solver_local(3, out_file);
	// 	test_prob.Solver_local_ext(3, out_file);
	// 	out_file << '\n';

	// }
	// out_file.close();
	
	// // Experiments for all test instances
	string file_path = "Dataset/MKIP/";
	string log_path = "Log/Log_MKIP";
	ofstream out_file;

	out_file.open("MKIP_results_m5.txt");
	int m1 = 1, m2 = 5;
	for (int n = 10; n <= 50; n+=10)
	{
		for (int k = 2; k < 7; k++)
		{
			for (int i = 0; i < 10; i++)
			{
				string file_name = "MKIP_n" + to_string(n) + "k" + to_string(k) + "j" + to_string(i) + "_m1_" + to_string(m1) + "_m2_" + to_string(m2);

				cout << file_name << endl;
				BMKP test_prob(file_name, file_path, log_path);

				//test_prob.WriteToMps(file_name);

				out_file << file_name << "\t";
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


	return 0;
}