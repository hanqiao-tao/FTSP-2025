// The authors thank Prof.Jordi Pereira for kindly providing us with the source code developed in Pereira, J. (2016). Procedures for the bin packing problem with precedence constraints. European Journal of Operational Research, 250, 794¨C806.
#include <stdio.h>
#include <iostream>
#include <math.h>
#include <string.h>
#include <vector>
#include <map>
#include <numeric>
#include <time.h>
#include <map>
#include <queue>
#include <set>
#include <stack>
#include <array>
#include <ilcplex/ilocplex.h>

using namespace std;

#ifndef CLK_TCK
#define CLK_TCK	CLOCKS_PER_SEC
#endif

#define  MALLOC(x,n,type) do									   \
{																   \
	if ((x = (type *) malloc( (n) * sizeof(type))) == NULL)		   \
	{															   \
	    fprintf(stderr,"out of memory\n");                         \
        fprintf(stderr,"x %d type\n",n);                           \
	    exit(1);                                                   \
	}															   \
} while (0)

typedef IloArray<IloNumArray> NumMatrix;
typedef IloArray<IloNumVarArray> NumVarMatrix;
#define NUM_JOB  1001  // 101
#define MAX_N  NUM_JOB + 1
#define BIG_INT 2147483647

// version: 0  FTSP solve;  1  sensitive analysis for parameter r;  2  sensitive analysis for parameter Ck 
extern int version;
extern unsigned int seed;							  
extern double intTolerance;

extern int nj;
extern int due_date;								  
extern int num_config;								 
extern vector<int> C_k;								
extern vector<int> C_k_sum;							
extern vector<int> Job_Time;					
extern vector<vector<bool>> config_matr;			 
extern vector<vector<bool>> adj_matr;               
extern vector<vector<bool>> trans_matr;            
extern vector<pair<int, int>> pre_rela;              
extern vector<vector<int>> direct_sucs;              
extern vector<vector<int>> direct_preds;             
extern vector<vector<int>> all_sucs;               
extern vector<vector<int>> all_preds;                 
extern vector<vector<int>> forbidden_set_2;			  
extern vector<vector<int>> forbidden_set_3; 		 
extern bool forbidden_set_matr_2[NUM_JOB][NUM_JOB];
extern bool forbidden_set_matr_3[NUM_JOB][NUM_JOB][NUM_JOB];
extern vector<int> configuration_sort_list;			

#define QS_COMPARE(a,b) ((b)-(a))
#define kMax 21  
extern int LB;
extern int LB_machine;
extern int LB_1;
extern int LB_DW;
extern double t_machine;
extern double t_lb_1;
extern double t_dw;
extern double t_pricing;
extern double t_el;
extern double t_H;

extern double* heads_c;
extern double* tails_c;
extern double* p;
extern int E[NUM_JOB];
extern int L[NUM_JOB];
extern vector<int> shortest_path;					

extern vector<vector<int>> H_solution;
extern double t_BDP;
extern int optimalBDP;								
extern int if_generate;								
extern int UB_BDP;
extern vector<vector<int>> BDP_solution;

extern double t_ILS;
extern vector<vector<int>> solution;
extern vector<vector<int>> best_solution;

struct Data
{
	double weight;
	Data* previous_ptr;
	char degrees[NUM_JOB];
};

struct compare_pq
{
	bool operator()(const Data& l, const Data& r)  const
	{
		return l.weight < r.weight;
	}
};

struct compare_set
{
	bool operator()(const Data& l, const Data& r)  const
	{
		return (memcmp(l.degrees, r.degrees, MAX_N) < 0);
	}
};

bool instance_reader_SALBP_otto(string filename);
bool reader_sort(string filename);
void FTSP_free();
void WriteLogFile(const string& szString, const string& filename);
void WriteLogFile_byEntry(const string& szString, const string& filename);
void Writeoutlog(string file_name);

vector<int> first_fit_BPP(vector<int> item_list);
void heuristic_H();

void SALBP_trans_closure_warshall();
void SALBP_matrix_to_graph_vec();
void compute_C();
void configuration_order();
int FTSP_solve();
int sensitive_r();
int sensitive_Ck();

void initMemoryLB4();
void freeMemoryLB4();
int get_pos(int k, int eps_pos, int task);
void Sort(int list[], double* pesos, int beg, int end, int posicionInicial);
void compute_data_LB4();
int compute_head(int task);
int compute_tail(int task);
void determine_Heads_or_Tails();
bool flow(int k_value, int task_value, int LB_try);
int determine_flow(int LB);
void compute_machine_lb();
void compute_dff_lbmodel(int kk, int eps_pos);
void compute_dff_lb();
double bound_pricing(int rest_time, int position);
int knapsack_branch(int position, int num_select);
void compute_forbidden_set();
void pricing_solve();
void compute_LB_DW();

bool dfs(int u);
bool hungarian();
void init_bounds();
void free_bounds();
int bound4(int ne, char* deg);
void initialize_bdp();
bool saveBestSolution();
void test_solution_complete();
int task_solitary(int n_eligible);
void gen_loads_bdp(int depth, int remaining_time, int start, int n_eligible, vector<int> config_fea);
void free_bdp();
void bdp();

void VNS_test_back();
void VNS_test_foward();
void VNS_flight();
void Pert();
void ILS_solve();
