/***********************************************************************************************************************
Xueyu Shi
University of Pittsburgh, Industrial Engineering
April 2020


Experinments for the bilevel minimum spanning tree problem
************************************************************************************************************************/
#pragma warning(disable: 4996)

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

// 3D variable array
typedef IloArray<IloNumVarArray> NumVarMatrix; 

class Graph
{
public:
	Graph()
		:num_edges_(0), num_nodes_(0) {}
	Graph(string filename, string file_path, string log_path){
		log_path_ = log_path;
		filename_ = filename;
		if (!ReadFile(filename, file_path))
			return;

		
		/* 
		*   Find the follower's controable edges 
		*/
		for (size_t i = 0; i < num_edges_; i++)
		{
			follower_edges_.push_back(i);
		}

		SortF(follower_edges_); // sort as the increasing order based on the optimality conditions
		vector<int> emptyset;

		vector<int> tmp_edges = Kruskal(emptyset); // Get E_f

		follower_edges_ = tmp_edges;

		fix_follower_edges_.resize(follower_edges_.size(), 0);
		
		ReduceLeader();

		int iter = 0;
		do
		{
			leader_edges_ = pre_leader_edges_;
			// Obtain the leader_nodes
			leader_nodes_.resize(num_nodes_);
			for (size_t i = 0; i < pre_leader_edges_.size(); i++)
			{

				leader_nodes_[vedge_[pre_leader_edges_[i]][0]] = 1;
				leader_nodes_[vedge_[pre_leader_edges_[i]][1]] = 1;
			}

			FixFollower();

			ReduceLeader();
			iter++;
		} while (leader_edges_.size() != pre_leader_edges_.size());

		
	}
	~Graph(){
		Clear();
	}

	void Clear(){
		num_edges_ = 0;
		num_nodes_ = 0;

		for (size_t i = 0; i < vedge_.size(); i++)
		{
			vedge_[i].clear();
		}
		vedge_.clear();

		for (size_t i = 0; i < nedge_.size(); i++)
		{
			for (size_t j = 0; j < nedge_[i].size(); j++)
			{
				nedge_[i][j].clear();
			}
			nedge_[i].clear();
		}
		nedge_.clear();

		c_.clear();
		d_.clear();
		leader_edges_.clear();
		follower_edges_.clear();

		leader_nodes_.clear();

		pre_leader_edges_.clear();
		fix_follower_edges_.clear();
	}

	bool ReadFile(string filename, string file_path){
		string line;
		ifstream in_file(file_path + "/" + filename + ".txt");

		Clear();

		if (in_file.is_open()) {
			int line_num = 0;
			while (getline(in_file, line)) {

				//cout << line << endl;
				

				vector<int> edge;
				vector<int> nedge(2);
				istringstream split(line);
				if (line_num == 0) {
					split >> num_nodes_ >> num_edges_;
					nedge_.resize(num_nodes_);
				}
				else if (line_num <= num_edges_)
				{
					for (string each; getline(split, each, '\t'); edge.push_back(stoi(each)));
					//				split >> edge[0] >> edge[1] >> edge[2] >> edge[3];
					c_.push_back(edge[2]);
					d_.push_back(edge[3]);
					edge[0] -= 1;
					edge[1] -= 1;

					vedge_.push_back(edge);

					nedge[0] = edge[0];
					nedge[1] = vedge_.size() - 1;
					nedge_[edge[0]].push_back(nedge);
					nedge[0] = edge[1];
					nedge_[edge[1]].push_back(nedge);
				}
				else
				{
					vector<string> split_line((istream_iterator<string>(split)), istream_iterator<string>());
					if (split_line.size() == 0) continue;
					for (size_t i = 1; i < split_line.size(); i++) leader_edges_.push_back(stoi(split_line[i]));
					//for (string each; getline(split, each, '\t'); leader_edges_.push_back(stoi(each)));
					for (size_t i = 0; i < leader_edges_.size(); i++)
					{
						leader_edges_[i] -= 1;
					}
				}
				line_num++;
			}


		}
		else {
			cout << "Can not open file " << filename << endl;
			return 0;
		}

		

		in_file.close();
		return 1;
	}

	


	bool Solver_EBP(ofstream& out_result){
		IloEnv env;	
		clock_t elpased_time = clock();
		try {
			IloModel model(env);

			// creat model global vairables
			IloIntVarArray xVars(env, 0); // leader's decision variable x
			IloIntVarArray yVars(env, 0); // follower's decision variable y

			for (size_t i = 0; i < pre_leader_edges_.size(); i++)
			{
				xVars.add(IloIntVar(env, 0, 1, create_name("x", i).c_str()));
			}

			for (size_t i = 0; i < follower_edges_.size(); i++)
			{
				yVars.add(IloIntVar(env, 0, 1, create_name("y", i).c_str()));
			}

			model.add(xVars);
			model.add(yVars);

			// define objective function
			IloExpr objExpr(env, 0);
			for (size_t i = 0; i < pre_leader_edges_.size(); i++)
			{
				objExpr += c_[pre_leader_edges_[i]] * xVars[i];
			}
			for (size_t i = 0; i < follower_edges_.size(); i++)
			{
				objExpr += c_[follower_edges_[i]] * yVars[i];
			}

			IloObjective obj(env, objExpr, IloObjective::Minimize, "objective function");
			model.add(obj);
			objExpr.end();

			/***************************************************
			* 		Constraints
	 		****************************************************/
			
			// constraint for interdiction
			for (size_t i = 0; i < pre_leader_edges_.size(); i++)
			{
				for (size_t j=0; j < follower_edges_.size(); j++)
				{
					if (pre_leader_edges_[i] == follower_edges_[j]) {
						model.add(xVars[i] + yVars[j] <= 1);
						break;
					}
				}
			}


			/* 
			* constraints to construct a spanning tree
			*/

			model.add(IloSum(xVars) + IloSum(yVars) == num_nodes_ - 1);
			
			// compute incidence matrix
			vector<vector<int>	> incident_matrix_leader(num_nodes_);

			for (size_t i = 0; i < pre_leader_edges_.size(); i++)
			{
				int node_1 = vedge_[pre_leader_edges_[i]][0];
				int node_2 = vedge_[pre_leader_edges_[i]][1]; // endpoints of the edge

				incident_matrix_leader[node_1].push_back(i);
				incident_matrix_leader[node_2].push_back(i);
			}

			vector<vector<int> > incident_matrix_follower(num_nodes_); // edges incident to the vertex, called incident matrix for the follower
			for (size_t i = 0; i < follower_edges_.size(); ++i)
			{
				int node_1 = vedge_[follower_edges_[i]][0];
				int node_2 = vedge_[follower_edges_[i]][1]; // endpoints of the edge

				incident_matrix_follower[node_1].push_back(i);
				incident_matrix_follower[node_2].push_back(i);

			}

			for (int i = 1; i < num_nodes_; ++i)
			{
				IloNumVarArray primal_xVars_0(env, 0);
				IloNumVarArray primal_yVars_0(env, 0);
				IloNumVarArray primal_xVars_1(env, 0); // note that it is a flow problem at undirected graph, hence we need to change edges to arcs with two directions. 
				IloNumVarArray primal_yVars_1(env, 0);

				for (size_t j = 0; j < pre_leader_edges_.size(); j++)
				{
					primal_xVars_0.add(IloNumVar(env, 0, +IloInfinity, create_name("f_l", i, j, 0).c_str()));
					primal_xVars_1.add(IloNumVar(env, 0, +IloInfinity, create_name("f_l", i, j, 1).c_str()));
				}
				for (size_t j = 0; j < follower_edges_.size(); j++)
				{
					primal_yVars_0.add(IloNumVar(env, 0, +IloInfinity, create_name("f_f", i, j, 0).c_str()));
					primal_yVars_1.add(IloNumVar(env, 0, +IloInfinity, create_name("f_f", i, j, 1).c_str()));
				}

				IloRangeArray flowCons(env, 0);
				for (size_t j = 0; j < pre_leader_edges_.size(); j++) {
					flowCons.add(primal_xVars_0[j] + primal_xVars_1[j] - xVars[j] <= 0);
				}
				for (size_t j = 0; j < follower_edges_.size(); j++)
				{
					flowCons.add(primal_yVars_0[j] + primal_yVars_1[j] - yVars[j] <= 0);
				}



				// constraints for flow constraint Af = 1, 0, -1
				for (int j = 0; j < num_nodes_; j++)
				{
					IloExpr flowExpr(env, 0);
					for (size_t k = 0; k < incident_matrix_leader[j].size(); k++)
					{
						if (j == vedge_[pre_leader_edges_[incident_matrix_leader[j][k]]][0]) // primal_0 direction from node 1 to node 2; primal_2 from node 2 to node 1
						{
							flowExpr += primal_xVars_0[incident_matrix_leader[j][k]] - primal_xVars_1[incident_matrix_leader[j][k]];
						}
						else
						{
							flowExpr += primal_xVars_1[incident_matrix_leader[j][k]] - primal_xVars_0[incident_matrix_leader[j][k]];
						}

					}

					for (size_t k = 0; k < incident_matrix_follower[j].size(); k++)
					{
						if (j == vedge_[follower_edges_[incident_matrix_follower[j][k]]][0])
							flowExpr += primal_yVars_0[incident_matrix_follower[j][k]] - primal_yVars_1[incident_matrix_follower[j][k]];
						else
						{
							flowExpr += primal_yVars_1[incident_matrix_follower[j][k]] - primal_yVars_0[incident_matrix_follower[j][k]];
						}
					}

					if (j == 0)
					{
						flowExpr -= 1;
					}
					if (j == i)
					{
						flowExpr += 1;
					}

					flowCons.add(flowExpr == 0);
					flowExpr.end();
				}

				model.add(flowCons);
				flowCons.end();
				primal_xVars_0.end();
				primal_xVars_1.end();
				primal_yVars_0.end();
				primal_yVars_1.end();
			}

			// constraints type 4
			int bigM = num_nodes_;
			vector<int> current_nodes = leader_nodes_;
			vector<vector<int> > incident_mat_l; // edges incident to the vertex, called incident matrix for the leader
			incident_mat_l.resize(num_nodes_);
			for (size_t i = 0; i < pre_leader_edges_.size(); i++)
			{
				int node_1 = vedge_[pre_leader_edges_[i]][0];
				int node_2 = vedge_[pre_leader_edges_[i]][1]; // endpoints of the edge

				incident_mat_l[node_1].push_back(i);
				incident_mat_l[node_2].push_back(i);
			}
			vector<vector<int> > incident_mat_f; // edges incident to the vertex, called incident matrix for the follower
			incident_mat_f.resize(num_nodes_);

			for (size_t i = 0; i < follower_edges_.size(); i++)
			{
				int node_1 = vedge_[follower_edges_[i]][0];
				int node_2 = vedge_[follower_edges_[i]][1]; // endpoints of the edge
				current_nodes[node_1] = 1;
				current_nodes[node_2] = 1; // current vertices

				if (fix_follower_edges_[i] == 1)
				{
					model.add(yVars[i] == 1);
					incident_mat_f[node_1].push_back(i);
					incident_mat_f[node_2].push_back(i);
					continue;
				}

				// constraints for dual variables
				IloNumVarArray piVars(env, 0);
				for (size_t j = 0; j < num_nodes_; j++)
				{
					piVars.add(IloNumVar(env, -IloInfinity, +IloInfinity, create_name("pi_f", i, j).c_str()));
				}
				model.add(piVars);

				IloRangeArray dualCons(env, 0);
				for (size_t j = 0; j < num_nodes_; j++)
				{
					if (current_nodes[j] == 0)
						dualCons.add(piVars[j] == 0);
				}

				// for node 1 and node 2
				dualCons.add(piVars[node_2] == 0);
				dualCons.add(piVars[node_1] - num_nodes_ * yVars[i] >= 0);
				dualCons.add(piVars[node_1] - yVars[i] <= num_nodes_ - 1);

				// for other nodes
				// edges in E_\ell
				for (size_t j = 0; j < pre_leader_edges_.size(); j++)
				{
					int tmp_node_1 = vedge_[pre_leader_edges_[j]][0];
					int tmp_node_2 = vedge_[pre_leader_edges_[j]][1];

					dualCons.add(piVars[tmp_node_1] - piVars[tmp_node_2] + bigM * xVars[j] <= bigM + 1);
					dualCons.add(piVars[tmp_node_2] - piVars[tmp_node_1] + bigM * xVars[j] <= bigM + 1);
				}
				// edges in E_f
				for (size_t j = 0; j < i; j++)
				{
					int tmp_node_1 = vedge_[follower_edges_[j]][0];
					int tmp_node_2 = vedge_[follower_edges_[j]][1];

					dualCons.add(piVars[tmp_node_1] - piVars[tmp_node_2] + bigM * yVars[j] <= bigM + 1);
					dualCons.add(piVars[tmp_node_2] - piVars[tmp_node_1] + bigM * yVars[j] <= bigM + 1);
				}

				model.add(dualCons);
				dualCons.end();
				

				// constraints for primal variables
				IloNumVarArray primal_xVars_0(env, 0);
				IloNumVarArray primal_yVars_0(env, 0);
				IloNumVarArray primal_xVars_1(env, 0); // note that it is a shortest path problem at undirected graph, hence we need to change edges to arcs with two directions. 
				IloNumVarArray primal_yVars_1(env, 0);
				IloNumVar	   primal_add(env, 0, +IloInfinity, create_name("x_a", i).c_str());

				for (size_t j = 0; j < pre_leader_edges_.size(); j++)
				{
					primal_xVars_0.add(IloNumVar(env, 0, +IloInfinity, create_name("x_f", i, j, 0).c_str()));
					primal_xVars_1.add(IloNumVar(env, 0, +IloInfinity, create_name("x_f", i, j, 1).c_str()));
				}
				for (size_t j = 0; j < i; j++)
				{
					primal_yVars_0.add(IloNumVar(env, 0, +IloInfinity, create_name("y_f", i, j, 0).c_str()));
					primal_yVars_1.add(IloNumVar(env, 0, +IloInfinity, create_name("y_f", i, j, 1).c_str()));
				}

				IloRangeArray primalCons(env, 0);
				for (size_t j = 0; j < pre_leader_edges_.size(); j++) {
					primalCons.add(primal_xVars_0[j] + primal_xVars_1[j] - xVars[j] <= 0);
				}
				for (size_t j = 0; j < i; j++)
				{
					primalCons.add(primal_yVars_0[j] + primal_yVars_1[j] - yVars[j] <= 0);
				}

				// constraints for shortest path Ax = 1, 0, -1
				for (size_t j = 0; j < num_nodes_; j++)
				{
					IloExpr flowCons(env, 0);
					for (size_t k = 0; k < incident_mat_l[j].size(); k++)
					{
						if (j == vedge_[pre_leader_edges_[incident_mat_l[j][k]]][0]) // primal_0 direction from node 1 to node 2; primal_2 from node 2 to node 1
						{
							flowCons += primal_xVars_0[incident_mat_l[j][k]] - primal_xVars_1[incident_mat_l[j][k]];
						}
						else
						{
							flowCons += primal_xVars_1[incident_mat_l[j][k]] - primal_xVars_0[incident_mat_l[j][k]];
						}

					}

					for (size_t k = 0; k < incident_mat_f[j].size(); k++)
					{
						if (j == vedge_[follower_edges_[incident_mat_f[j][k]]][0])
							flowCons += primal_yVars_0[incident_mat_f[j][k]] - primal_yVars_1[incident_mat_f[j][k]];
						else
						{
							flowCons += primal_yVars_1[incident_mat_f[j][k]] - primal_yVars_0[incident_mat_f[j][k]];
						}
					}

					if (j == node_1)
					{
						flowCons += 1 - primal_add;
					}
					if (j == node_2)
					{
						flowCons += primal_add - 1;
					}

					primalCons.add(flowCons == 0);
				}
				model.add(primalCons);
				primalCons.end();

				// constraints for strong duality
				IloExpr opt_expr(env, 0);
				for (size_t j = 0; j < pre_leader_edges_.size(); j++)
				{
					opt_expr += primal_xVars_0[j] + primal_xVars_1[j];
				}

				for (size_t j = 0; j < i; j++)
				{
					opt_expr += primal_yVars_0[j] + primal_yVars_1[j];
				}
				opt_expr += num_nodes_ * primal_add - piVars[node_1];
				model.add(opt_expr == 0);
				opt_expr.end();

				primal_xVars_0.end();
				primal_xVars_1.end();
				primal_yVars_0.end();
				primal_yVars_1.end();
				piVars.end();


				// update incident matrix
				incident_mat_f[node_1].push_back(i);
				incident_mat_f[node_2].push_back(i);
				//printMat(incident_mat_f);
			}

			// solve 
			IloCplex cplex(env);
			cplex.extract(model);
			cplex.setParam(IloCplex::TiLim, 3600);
			cplex.setParam(IloCplex::Param::Threads, 1);
			string cplex_log = log_path_ + "/cplex_" + filename_ + "_ext.log";
			ofstream out_file(cplex_log, std::ios_base::out);
			cplex.setOut(out_file);

			//cplex.exportModel("MIPM_E.lp");

			
			//cplex.exportModel("model.lp");
			cplex.solve();
			elpased_time = clock() - elpased_time;

			double mip_obj = cplex.getObjValue();

			
			// get values
			/*IloNumArray xOpt(env);
			cplex.getValues(xOpt, xVars);
			IloNumArray yOpt(env);
			cplex.getValues(yOpt, yVars);

			env.out() << "Leader's optimal decision: " << xOpt << endl;
			cout << "Follower's edges: " << endl;
			printVec(&follower_edges_);
			env.out() << "Follower's optimal decision: " << yOpt << endl;*/

			// lp relaxation
			for (size_t i = 0; i < pre_leader_edges_.size(); i++)
			{
				model.add(IloConversion(env, xVars[i], ILOFLOAT));
			}

			for (size_t i = 0; i < follower_edges_.size(); i++)
			{
				model.add(IloConversion(env, yVars[i], ILOFLOAT));
			}
			cplex.extract(model);

			//cplex.exportModel("LPMIP.lp");
			cplex.solve();

			double lp_obj = cplex.getObjValue();
			double gap = (mip_obj - lp_obj) / mip_obj * 100;

			out_result << lp_obj << "\t" << mip_obj << "\t" << gap << "\t" << (float)elpased_time / CLOCKS_PER_SEC << '\t';

			cout <<  lp_obj << "\t" << mip_obj << "\t" << gap << "\t" << (float)elpased_time / CLOCKS_PER_SEC << endl;
		
			model.end();
			env.end();

			return 1;


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

	bool Solver_BP(ofstream& out_result){
		IloEnv env;	
		clock_t elpased_time = clock();
		try {
			IloModel model(env);

			// creat model global vairables
			IloIntVarArray xVars(env, 0); // leader's decision variable x
			IloIntVarArray yVars(env, 0); // follower's decision variable y

			for (size_t i = 0; i < pre_leader_edges_.size(); i++)
			{
				xVars.add(IloIntVar(env, 0, 1, create_name("x", i).c_str()));
			}

			for (size_t i = 0; i < follower_edges_.size(); i++)
			{
				yVars.add(IloIntVar(env, 0, 1, create_name("y", i).c_str()));
			}

			model.add(xVars);
			model.add(yVars);

			// define objective function
			IloExpr objExpr(env, 0);
			for (size_t i = 0; i < pre_leader_edges_.size(); i++)
			{
				objExpr += c_[pre_leader_edges_[i]] * xVars[i];
			}
			for (size_t i = 0; i < follower_edges_.size(); i++)
			{
				objExpr += c_[follower_edges_[i]] * yVars[i];
			}

			IloObjective obj(env, objExpr, IloObjective::Minimize, "objective function");
			model.add(obj);
			objExpr.end();

			/***************************************************
			* 		Constraints
	 		****************************************************/
			
			// constraint for interdiction
			for (size_t i = 0; i < pre_leader_edges_.size(); i++)
			{
				for (size_t j=0; j < follower_edges_.size(); j++)
				{
					if (pre_leader_edges_[i] == follower_edges_[j]) {
						model.add(xVars[i] + yVars[j] <= 1);
						break;
					}
				}
			}

			/* 
			* constraints to construct a spanning tree
			*/

			model.add(IloSum(xVars) + IloSum(yVars) == num_nodes_ - 1);
			
			// compute incidence matrix
			vector<vector<int>	> incident_matrix_leader(num_nodes_);

			for (size_t i = 0; i < pre_leader_edges_.size(); i++)
			{
				int node_1 = vedge_[pre_leader_edges_[i]][0];
				int node_2 = vedge_[pre_leader_edges_[i]][1]; // endpoints of the edge

				incident_matrix_leader[node_1].push_back(i);
				incident_matrix_leader[node_2].push_back(i);
			}

			vector<vector<int> > incident_matrix_follower(num_nodes_); // edges incident to the vertex, called incident matrix for the follower
			for (size_t i = 0; i < follower_edges_.size(); ++i)
			{
				int node_1 = vedge_[follower_edges_[i]][0];
				int node_2 = vedge_[follower_edges_[i]][1]; // endpoints of the edge

				incident_matrix_follower[node_1].push_back(i);
				incident_matrix_follower[node_2].push_back(i);

			}

			for (int i = 1; i < num_nodes_; ++i)
			{
				IloNumVarArray primal_xVars_0(env, 0);
				IloNumVarArray primal_yVars_0(env, 0);
				IloNumVarArray primal_xVars_1(env, 0); // note that it is a flow problem at undirected graph, hence we need to change edges to arcs with two directions. 
				IloNumVarArray primal_yVars_1(env, 0);

				for (size_t j = 0; j < pre_leader_edges_.size(); j++)
				{
					primal_xVars_0.add(IloNumVar(env, 0, +IloInfinity, create_name("f_l", i, j, 0).c_str()));
					primal_xVars_1.add(IloNumVar(env, 0, +IloInfinity, create_name("f_l", i, j, 1).c_str()));
				}
				for (size_t j = 0; j < follower_edges_.size(); j++)
				{
					primal_yVars_0.add(IloNumVar(env, 0, +IloInfinity, create_name("f_f", i, j, 0).c_str()));
					primal_yVars_1.add(IloNumVar(env, 0, +IloInfinity, create_name("f_f", i, j, 1).c_str()));
				}

				IloRangeArray flowCons(env, 0);
				for (size_t j = 0; j < pre_leader_edges_.size(); j++) {
					flowCons.add(primal_xVars_0[j] + primal_xVars_1[j] - xVars[j] <= 0);
				}
				for (size_t j = 0; j < follower_edges_.size(); j++)
				{
					flowCons.add(primal_yVars_0[j] + primal_yVars_1[j] - yVars[j] <= 0);
				}



				// constraints for flow constraint Af = 1, 0, -1
				for (int j = 0; j < num_nodes_; j++)
				{
					IloExpr flowExpr(env, 0);
					for (size_t k = 0; k < incident_matrix_leader[j].size(); k++)
					{
						if (j == vedge_[pre_leader_edges_[incident_matrix_leader[j][k]]][0]) // primal_0 direction from node 1 to node 2; primal_2 from node 2 to node 1
						{
							flowExpr += primal_xVars_0[incident_matrix_leader[j][k]] - primal_xVars_1[incident_matrix_leader[j][k]];
						}
						else
						{
							flowExpr += primal_xVars_1[incident_matrix_leader[j][k]] - primal_xVars_0[incident_matrix_leader[j][k]];
						}

					}

					for (size_t k = 0; k < incident_matrix_follower[j].size(); k++)
					{
						if (j == vedge_[follower_edges_[incident_matrix_follower[j][k]]][0])
							flowExpr += primal_yVars_0[incident_matrix_follower[j][k]] - primal_yVars_1[incident_matrix_follower[j][k]];
						else
						{
							flowExpr += primal_yVars_1[incident_matrix_follower[j][k]] - primal_yVars_0[incident_matrix_follower[j][k]];
						}
					}

					if (j == 0)
					{
						flowExpr -= 1;
					}
					if (j == i)
					{
						flowExpr += 1;
					}

					flowCons.add(flowExpr == 0);
					flowExpr.end();
				}

				model.add(flowCons);
				flowCons.end();
				primal_xVars_0.end();
				primal_xVars_1.end();
				primal_yVars_0.end();
				primal_yVars_1.end();
			}

			// constraints type 4
			int bigM = num_nodes_;
			for (int i = 0; i < follower_edges_.size(); i++)
			{
				if (fix_follower_edges_[i] == 1)
				{
					model.add(yVars[i] == 1);
					continue;
				}

				int node_1 = vedge_[follower_edges_[i]][0];
				int node_2 = vedge_[follower_edges_[i]][1];

				for (int j = i+1; j < follower_edges_.size(); j++)
				{
					IloNumVarArray primal_xVars_0(env, 0);
					IloNumVarArray primal_yVars_0(env, 0);
					IloNumVarArray primal_xVars_1(env, 0); // note that it is a flow problem at undirected graph, hence we need to change edges to arcs with two directions. 
					IloNumVarArray primal_yVars_1(env, 0);
					IloNumVar primal_add(env, 0, +IloInfinity, create_name("z_", i, j).c_str());

					for (size_t k = 0; k < pre_leader_edges_.size(); k++)
					{
						primal_xVars_0.add(IloNumVar(env, 0, +IloInfinity, create_name("z_l0_", k, i, j).c_str()));
						primal_xVars_1.add(IloNumVar(env, 0, +IloInfinity, create_name("z_l1_", k, i, j).c_str()));
					}
					for (size_t k = 0; k < follower_edges_.size(); k++)
					{
						primal_yVars_0.add(IloNumVar(env, 0, +IloInfinity, create_name("z_f0_", k, i, j).c_str()));
						primal_yVars_1.add(IloNumVar(env, 0, +IloInfinity, create_name("z_f1_", k, i, j).c_str()));
					}

					IloRangeArray flowCons(env, 0);
					for (size_t k = 0; k < pre_leader_edges_.size(); k++) {
						flowCons.add(primal_xVars_0[k] + primal_xVars_1[k] - xVars[k] <= 0);
					}
					for (size_t k = 0; k < follower_edges_.size(); k++)
					{
						if(k==i || k==j) {
							flowCons.add(primal_yVars_0[k] == 0);
							flowCons.add(primal_yVars_1[k] == 0);
							continue;
						}
						flowCons.add(primal_yVars_0[k] + primal_yVars_1[k] - yVars[k] <= 0);
					}



					// constraints for flow constraint Af = 1, 0, -1
					for (int k = 0; k < num_nodes_; k++)
					{
						IloExpr flowExpr(env, 0);
						for (size_t l = 0; l < incident_matrix_leader[k].size(); l++)
						{
							if (k == vedge_[pre_leader_edges_[incident_matrix_leader[k][l]]][0]) // primal_0 direction from node 1 to node 2; primal_2 from node 2 to node 1
							{
								flowExpr += primal_xVars_0[incident_matrix_leader[k][l]] - primal_xVars_1[incident_matrix_leader[k][l]];
							}
							else
							{
								flowExpr += primal_xVars_1[incident_matrix_leader[k][l]] - primal_xVars_0[incident_matrix_leader[k][l]];
							}

						}

						for (size_t l = 0; l < incident_matrix_follower[k].size(); l++)
						{
							if ((incident_matrix_follower[k][l] == i) || (incident_matrix_follower[k][l] == j))
							{
								continue;
							}

							if (k == vedge_[follower_edges_[incident_matrix_follower[k][l]]][0])
								flowExpr += primal_yVars_0[incident_matrix_follower[k][l]] - primal_yVars_1[incident_matrix_follower[k][l]];
							else
							{
								flowExpr += primal_yVars_1[incident_matrix_follower[k][l]] - primal_yVars_0[incident_matrix_follower[k][l]];
							}
						}

						if (k == node_1)
						{
							flowExpr += primal_add - 1;
						}
						if (k == node_2)
						{
							flowExpr += 1 - primal_add;
						}

						flowCons.add(flowExpr == 0);
						flowExpr.end();
					}

					model.add(flowCons);
					flowCons.end();

					model.add(IloSum(primal_xVars_0) + IloSum(primal_xVars_1) + IloSum(primal_yVars_0) + IloSum(primal_yVars_1) + num_nodes_*primal_add -yVars[i] + yVars[j] <= num_nodes_);
					primal_xVars_0.end();
					primal_xVars_1.end();
					primal_yVars_0.end();
					primal_yVars_1.end();
				}
			}

			// solve 
			IloCplex cplex(env);
			cplex.extract(model);
			cplex.setParam(IloCplex::TiLim, 3600);
			cplex.setParam(IloCplex::Param::Threads, 1);
			string cplex_log = log_path_ + "/cplex_" + filename_ + "_BP.log";
			ofstream out_file(cplex_log, std::ios_base::out);
			cplex.setOut(out_file);

			//cplex.exportModel("MIPM_E.lp");

			cplex.solve();
			elpased_time = clock() - elpased_time;

			double mip_obj = cplex.getObjValue();

			
			// get values
			/*IloNumArray xOpt(env);
			cplex.getValues(xOpt, xVars);
			IloNumArray yOpt(env);
			cplex.getValues(yOpt, yVars);

			env.out() << "Leader's optimal decision: " << xOpt << endl;
			cout << "Follower's edges: " << endl;
			printVec(&follower_edges_);
			env.out() << "Follower's optimal decision: " << yOpt << endl;*/

			// lp relaxation
			for (size_t i = 0; i < pre_leader_edges_.size(); i++)
			{
				model.add(IloConversion(env, xVars[i], ILOFLOAT));
			}

			for (size_t i = 0; i < follower_edges_.size(); i++)
			{
				model.add(IloConversion(env, yVars[i], ILOFLOAT));
			}
			cplex.extract(model);

			//cplex.exportModel("LPMIP.lp");
			cplex.solve();

			double lp_obj = cplex.getObjValue();
			double gap = (mip_obj - lp_obj) / mip_obj * 100;

			out_result << num_nodes_ << "\t" << num_edges_ << "\t" << lp_obj << "\t" << mip_obj << "\t" << gap << "\t" << (float)elpased_time / CLOCKS_PER_SEC << '\t';

			cout << num_nodes_ << "\t" << num_edges_ <<  "\t" << lp_obj << "\t" << mip_obj << "\t" << gap << "\t" << (float)elpased_time / CLOCKS_PER_SEC << '\t';

			
			model.end();
			env.end();
			

			return 1;


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
	// preprocessing
	void ReduceLeader(){
		int max_weight_follower = 0;
		for (size_t i = 0; i < follower_edges_.size(); i++)
		{
			if (fix_follower_edges_[i] == 1)
			{
				continue;
			}
			if (c_[follower_edges_[i]] > max_weight_follower)
				max_weight_follower = c_[follower_edges_[i]];
		}

		pre_leader_edges_.resize(0);
		for (size_t i = 0; i < leader_edges_.size(); i++)
		{
			if (c_[leader_edges_[i]] < max_weight_follower)
			{
				pre_leader_edges_.push_back(leader_edges_[i]);
			}
		}
	}
	void FixFollower(){
		// bulid the component system based on E_l
		vector<vector<int> > component_sets; // the connected components in the tree
		int num_components = num_nodes_;
		for (int i = 0; i < num_nodes_; i++) { // initialize the component by the individual nodes
			vector<int> node_set;
			node_set.push_back(i);
			component_sets.push_back(node_set);
		}
		vector<int> node_index; // indicate the component index of the node
		for (size_t i = 0; i < num_nodes_; i++)
		{
			node_index.push_back(i);

		}

		for (size_t i = 0; i < pre_leader_edges_.size(); i++)
		{
			int node_1 = vedge_[pre_leader_edges_[i]][0];
			int node_2 = vedge_[pre_leader_edges_[i]][1];

			if (node_index[node_1] == node_index[node_2])
				continue;

			int set_num_1 = node_index[node_1];
			int set_num_2 = node_index[node_2];
			for (size_t j = 0; j < component_sets[set_num_2].size(); j++)
			{
				node_index[component_sets[set_num_2][j]] = set_num_1;
			}
			component_sets[set_num_1].insert(component_sets[set_num_1].end(), component_sets[set_num_2].begin(), component_sets[set_num_2].end());
			component_sets[set_num_2].clear();

			num_components -= 1;
		}

		//

		fix_follower_edges_.clear();
		fix_follower_edges_.resize(follower_edges_.size());
		for (size_t i = 0; i < follower_edges_.size(); i++)
		{
			if (num_components == 1) {
				return;
			}
			// endpoints
			int node_1 = vedge_[follower_edges_[i]][0];
			int node_2 = vedge_[follower_edges_[i]][1];

			if (node_index[node_1] == node_index[node_2])
				continue;

			// fix edge
			fix_follower_edges_[i] = 1;

			// merge 
			int set_num_1 = node_index[node_1];
			int set_num_2 = node_index[node_2];
			for (size_t j = 0; j < component_sets[set_num_2].size(); j++)
			{
				node_index[component_sets[set_num_2][j]] = set_num_1;
			}
			component_sets[set_num_1].insert(component_sets[set_num_1].end(), component_sets[set_num_2].begin(), component_sets[set_num_2].end());
			component_sets[set_num_2].clear();

			num_components -= 1;
		}
	}

	// define the operator between nodes under the optimistic rule
	bool edge_order(int i, int j){
		return((d_[i] < d_[j]) || (d_[i] == d_[j]) && (c_[i] < c_[j]));
	}
	// Solve the reaction problem with resepect to the leader's decision by Kruskal algorithm
	vector<int> Kruskal(vector<int> leader_decision){
		vector<int> follower_decison; // the decision taken by the follower

		vector<vector<int> > component_sets; // the connected components in the tree
		for (int i = 0; i < num_nodes_; i++) { // initialize the component by the individual nodes
			vector<int> node_set;
			node_set.push_back(i);
			component_sets.push_back(node_set);
		}
		vector<int> node_index; // indicate the component index of the node
		for (size_t i = 0; i < num_nodes_; i++)
		{
			node_index.push_back(i);
		}


		// We first fix the leader decision into the components
		for (size_t i = 0; i < leader_decision.size(); i++)
		{
			// endpoints of the edge
			int node_1 = vedge_[leader_decision[i]][0];
			int node_2 = vedge_[leader_decision[i]][1];

			if (node_index[node_1] == node_index[node_2])
			{
				cout << "L is not feasible." << endl;
				return follower_decison;
			}

			// merge the components of node_1 and node_2
			int set_num_1 = node_index[node_1];
			int set_num_2 = node_index[node_2];
			for (size_t i = 0; i < component_sets[set_num_2].size(); i++)
			{
				node_index[component_sets[set_num_2][i]] = set_num_1;
			}
			component_sets[set_num_1].insert(component_sets[set_num_1].end(), component_sets[set_num_2].begin(), component_sets[set_num_2].end());
			component_sets[set_num_2].clear();
		}

		// Decides edges controlled by the follower
		for (size_t i = 0; i < follower_edges_.size(); i++)
		{
			// endpoints of the edge
			int node_1 = vedge_[follower_edges_[i]][0];
			int node_2 = vedge_[follower_edges_[i]][1];

			if (node_index[node_1] == node_index[node_2]) continue; // drop the edge


			follower_decison.push_back(follower_edges_[i]); // add the edge into the follower's decision

															// merge the components of node_1 and node_2
			int set_num_1 = node_index[node_1];
			int set_num_2 = node_index[node_2];
			for (size_t i = 0; i < component_sets[set_num_2].size(); i++)
			{
				node_index[component_sets[set_num_2][i]] = set_num_1;
			}
			component_sets[set_num_1].insert(component_sets[set_num_1].end(), component_sets[set_num_2].begin(), component_sets[set_num_2].end());
			component_sets[set_num_2].clear();

			// already a spanning tree

			if (follower_decison.size() + leader_decision.size() >= (num_nodes_ - 1)) {
				//cout << "break cond" << endl;
				break;
			}
		}

		return follower_decison;
	}

	// sort the edges over edge_oder() as the nondecreasing order
	void SortF(vector<int> &edges){
		int tmp;
		for (size_t i = 0; i < edges.size(); i++)
		{
			// find the lowest i th index
			for (size_t j = edges.size() - 1; j > i; j--)
			{
				if (edge_order(edges[j], edges[j - 1])) {
					tmp = edges[j - 1];
					edges[j - 1] = edges[j];
					edges[j] = tmp;
				}
			}

		}
	}

	string create_name(string letters, int number1, int number2 = -1, int number3 = -1){
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


	string filename_;
	string log_path_;
	int num_edges_;
	int num_nodes_;

	vector<vector<int> > vedge_; // vector of edges 
	vector<vector<vector<int> > > nedge_; // node-node edge matrix

	vector<int> c_; // edge cost of the leader
	vector<int> d_; // edge cost of the follower

	vector<int> leader_edges_; // edges controlled by the leader 
	vector<int> follower_edges_; // edges constrolled by the follower

	vector<int> leader_nodes_;

	vector<int> pre_leader_edges_; // edges controlled by the leader after preprocessing
	vector<int> fix_follower_edges_; // edges contained in the follower's reation solution whenever the leader's decision

};



int main() {

		
	// solve the single file
	/*ofstream out_file("results.txt");

	vector<double> result;
	Graph test_g("b01_1_0.3_8", "B/advanced_B", result);

	//test_g.WriteFile("output.txt");
	test_g.Solver(out_file);

	test_g.Solver_EBP(out_file);
	test_g.Solver_BP(out_file);*/


	ofstream out_file("BMST_results.txt");
	string file_path = "Dataset/BMST/";
	string log_path = "Log/Log_BMST";

	out_file << "\t INS. \t \t | \t \t BMST \t \t | \t \t EBMST \t \t " << endl;
	out_file << "File \t |N| \t |E| | \t LP \t MIP \t GAP \t Time | \t LP \t MIP \t GAP \t Time" << endl;

	for (int i = 1; i <= 18; ++i)
	{
		for (int j = 5; j <= 15; j+=5)
		{
			for (int k = 1; k <= 10; ++k)
			{
				string file_name;
				if (i < 10)
				{
					file_name = "b0" + to_string(i) + "_0_" + to_string(j) + "_" + to_string(k);
				}else{
					file_name = "b" + to_string(i) + "_0_" + to_string(j) + "_" + to_string(k);

				}

				cout << file_name << '\t';
				out_file << file_name << "\t";
				Graph test_graph(file_name, file_path, log_path);
				

				test_graph.Solver_BP(out_file);
				test_graph.Solver_EBP(out_file);
				
				out_file << '\n';
			}
		}
		
	}

	out_file.close();


	return 0;
}
