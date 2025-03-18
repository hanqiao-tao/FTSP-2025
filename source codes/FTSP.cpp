#include "FTSP.h"

unsigned int seed = 922;
double intTolerance = 0.0001;
int if_generate = 1;
int version = 0;

int nj;
int due_date;
int num_config;
vector<int> C_k;
vector<int> C_k_sum;
vector<int> Job_Time;
vector<vector<bool>> config_matr;
vector<vector<bool>> adj_matr;
vector<vector<bool>> trans_matr;
vector<pair<int, int>> pre_rela;
vector<vector<int>> direct_sucs;
vector<vector<int>> direct_preds;
vector<vector<int>> all_sucs;
vector<vector<int>> all_preds;
vector<int> positional_weight;
vector<vector<int>> forbidden_set_2;
vector<vector<int>> forbidden_set_3;
bool forbidden_set_matr_2[NUM_JOB][NUM_JOB];
bool forbidden_set_matr_3[NUM_JOB][NUM_JOB][NUM_JOB];
vector<int> configuration_sort_list;

int LB;
int LB_machine;
int LB_1;
int LB_DW;
double t_machine;
double t_lb_1;
double t_dw;
double t_pricing;
double t_el;
double t_H;

double* heads_c;
double* tails_c;
double* p;
int E[NUM_JOB];
int L[NUM_JOB];
vector<int> shortest_path;

vector<vector<int>> H_solution;

double t_BDP;
int optimalBDP;
int UB_BDP;
vector<vector<int>> BDP_solution;

double t_ILS;
vector<vector<int>> solution;
vector<vector<int>> best_solution;

#pragma region Data preparation
void SALBP_trans_closure_warshall()
{
    trans_matr = adj_matr;
    for (int k = 0; k <= nj; k++)
    {
        for (int i = 0; i <= nj; i++)
        {
            for (int j = 0; j <= nj; j++)
            {
                trans_matr[i][j] = trans_matr[i][j] || (trans_matr[i][k] && trans_matr[k][j]);
            }
        }
    }
}

void SALBP_matrix_to_graph_vec()
{
    pre_rela.reserve(nj * nj / 2);
    for (int i = 0; i < nj; i++)
    {
        for (int j = i + 1; j <= nj; j++)
        {
            if (adj_matr[i][j])
            {
                pre_rela.push_back(make_pair(i, j));
            }
        }
    }
    all_sucs.resize(nj + 1);
    all_preds.resize(nj + 1);
    direct_sucs.resize(nj + 1);
    direct_preds.resize(nj + 1);
    for (int i = 0; i <= nj; i++)
    {
        all_sucs[i].reserve(nj);
        all_preds[i].reserve(nj);
    }
    for (int i = 0; i <= nj; i++)
    {
        for (int j = i + 1; j <= nj; j++)
        {
            if (trans_matr[i][j])
                all_sucs[i].push_back(j);
            if (adj_matr[i][j])
                direct_sucs[i].push_back(j);
        }
        for (int j = 0; j < i; j++)
        {
            if (trans_matr[j][i])
                all_preds[i].push_back(j);
            if (adj_matr[j][i])
                direct_preds[i].push_back(j);
        }
    }
    shortest_path.resize(nj + 1, 0);
    int count = 1;
    while (count != 0)
    {
        count = 0;
        for (int i = 1; i <= nj; i++)
        {
            for (int j = 1; j <= nj; j++)
            {
                if ((i != j) && (trans_matr[j][i] == 1) && (shortest_path[i] < shortest_path[j]))
                {
                    shortest_path[i] = shortest_path[j];
                    count = 1;
                }
                if ((i != j) && (trans_matr[j][i] == 1) && (shortest_path[j] == shortest_path[i]))
                {
                    shortest_path[i] = shortest_path[j] + 1;
                    count = 1;
                }
            }
        }
    }
    positional_weight.resize(nj + 1);
    for (int i = 1; i <= nj; i++)
    {
        int sum = Job_Time[i];
        for (int j = i + 1; j <= nj; j++)
        {
            if (trans_matr[i][j])
            {
                sum = sum + Job_Time[j];
            }
        }
        positional_weight[i] = sum;
    }
}

void compute_C()
{
    C_k = vector<int>(num_config + 1);
    int L = 0;
    for (int i = 1; i <= nj; i++)
    {
        if (shortest_path[i] > L)
        {
            L = shortest_path[i];
        }
    }
    for (int k = 1; k <= num_config; k++)
    {
        int sum_Ik_t = 0;
        int nj_k = 0;
        for (int i = 1; i <= nj; i++)
        {
            if (config_matr[k][i])
            {
                sum_Ik_t += Job_Time[i];
                nj_k++;
            }
        }
        C_k[k] = L + 1 + (int)(ceil(2 * sum_Ik_t / due_date - intTolerance));
        C_k[k] = C_k[k] > nj_k ? nj_k : C_k[k];
    }
    vector<vector<vector<int>>>task_set_jk(num_config + 1);
    for (int k = 1; k <= num_config; k++)
    {
        task_set_jk[k].resize(L + 1);
        for (int j = 0; j <= L; j++)
        {
            task_set_jk[k][j].reserve(nj + 1);
            task_set_jk[k][j].push_back(0); 
        }
    }
    for (int i = 1; i <= nj; i++)
    {
        for (int k = 1; k <= num_config; k++)
        {
            if (config_matr[k][i])
            {
                int j = shortest_path[i];
                task_set_jk[k][j].push_back(i);
            }
        }
    }
    for (int k = 1; k <= num_config; k++)
    {
        int C_k_new = 0;
        for (int j = 0; j <= L; j++)
        {
            C_k_new += first_fit_BPP(task_set_jk[k][j])[0];
        }
        C_k[k] = C_k[k] > C_k_new ? C_k_new : C_k[k];
    }
    C_k_sum = vector<int>(num_config + 1, 0);
    for (int k = 1; k <= num_config; k++)
    {
        for (int kk = 1; kk <= k; kk++)
        {
            C_k_sum[k] += C_k[kk];
        }
    }

}

void configuration_order()
{
    configuration_sort_list.resize(num_config + 1);
    for (int k = 1; k <= num_config; k++) configuration_sort_list[k] = k;
    vector<int> config_positional_weight(num_config + 1);
    for (int i = 1; i <= nj; i++)
    {
        for (int k = 1; k <= num_config; k++)
        {
            if (config_matr[k][i]) config_positional_weight[k] += positional_weight[i];
        }
    }
    int change = 1;
    while (change > 0)
    {
        change = -1;
        for (int k = 1; k <= num_config; k++)
        {
            if (config_positional_weight[configuration_sort_list[k]] < config_positional_weight[configuration_sort_list[k]])
            {
                change = configuration_sort_list[k];
                configuration_sort_list[k] = configuration_sort_list[k + 1];
                configuration_sort_list[k + 1] = change;
            }
        }
    }
}

void FTSP_free()
{
    vector<int>().swap(Job_Time);
    vector<vector<bool>>().swap(adj_matr);
    vector<vector<bool>>().swap(trans_matr);
    vector<pair<int, int>>().swap(pre_rela);
    vector<vector<int>>().swap(direct_sucs);
    vector<vector<int>>().swap(direct_preds);
    vector<vector<int>>().swap(all_sucs);
    vector<vector<int>>().swap(all_preds);
    vector<int>().swap(shortest_path);
    vector<vector<bool>>().swap(config_matr);
    vector<int>().swap(positional_weight);
    vector<vector<int>>().swap(H_solution);
    vector<vector<int>>().swap(BDP_solution);
    vector<vector<int>>().swap(solution);
    vector<vector<int>>().swap(best_solution);
    vector<vector<int>>().swap(forbidden_set_2);
    vector<vector<int>>().swap(forbidden_set_3);
    vector<int>().swap(configuration_sort_list);

    for (int i = 1; i <= nj; i++)
    {
        memset(forbidden_set_matr_2[i], 0, sizeof(forbidden_set_matr_2[i]));
    }
    for (int i = 1; i <= nj; i++)
    {
        for (int j = 1; j <= nj; j++)
        {
            memset(forbidden_set_matr_3[i][j], 0, sizeof(forbidden_set_matr_3[i][j]));
        }
    }
    memset(E, 0, sizeof(E));
    memset(L, 0, sizeof(L));
    t_pricing = 0;
    if_generate = 1;
    t_BDP = 0.0;
    t_dw = 0.0;
    t_lb_1 = 0.0;
    t_ILS = 0.0;
    t_machine = 0.0;
    t_pricing = 0.0;
    t_el = 0.0;
    t_H = 0.0;
}
#pragma endregion

int FTSP_solve()
{
    string out_name = ".\\result.txt";
    WriteLogFile("", out_name);
    WriteLogFile_byEntry("LB", out_name);
    WriteLogFile_byEntry("LB_DFF", out_name);
    WriteLogFile_byEntry("LB_machine", out_name);
    WriteLogFile_byEntry("LB_DW", out_name);
    WriteLogFile_byEntry("UB", out_name);
    WriteLogFile_byEntry("initial solution", out_name);
    WriteLogFile_byEntry("if_BDP_optimal", out_name);
    WriteLogFile_byEntry("t_dff", out_name);
    WriteLogFile_byEntry("t_machine", out_name);
    WriteLogFile_byEntry("t_dw", out_name);
    WriteLogFile_byEntry("t_el", out_name);
    WriteLogFile_byEntry("t_H", out_name);
    WriteLogFile_byEntry("t_BDP", out_name);
    WriteLogFile_byEntry("t_ILS", out_name);
    WriteLogFile_byEntry("t_pricing", out_name);
    vector<int> ins_dd;
    for (int i = 1; i <= 525; i++)
    {
        ins_dd.push_back(i);
    }
    for (int ins_c = 0; ins_c < ins_dd.size(); ins_c++)
    {
        bool read_ins = instance_reader_SALBP_otto(".\\Otto dataset\\medium data set_n=50\\instance_n=50_" + to_string(ins_dd[ins_c]) + ".alb");
        cout << "instance_n = " << nj << "_" << ins_dd[ins_c] << ".alb " << endl;

        SALBP_trans_closure_warshall();
        SALBP_matrix_to_graph_vec();

        compute_machine_lb();
        for (num_config = 2; num_config <= 2; num_config++)
        {
            string filename = ".\\FTSP dataset\\S_" + to_string(nj) + "_" + to_string(num_config) + "_" + to_string(1) + ".txt";
            reader_sort(filename);
            compute_C();
            configuration_order();
            heuristic_H();
            compute_forbidden_set();
            bdp();
            if (solution[0][0] > H_solution[0][0])
            {
                solution = H_solution;
            }

            if (if_generate)
            {
                best_solution = solution;
                compute_dff_lb();
                compute_LB_DW();
                LB = LB_machine > LB_1 ? LB_machine : LB_1;
                LB = LB_DW > LB ? LB_DW : LB;
                if (LB < solution[0][0] && optimalBDP != 1) ILS_solve();
                else LB = solution[0][0];
                cout << best_solution[0][0] << "\t" << LB << endl;
            }
            else
            {
                cout << "no feasible solution found!" << endl;
            }
            Writeoutlog(out_name);
            vector<vector<int>>().swap(forbidden_set_2);
            vector<vector<int>>().swap(forbidden_set_3);
            for (int i = 1; i <= nj; i++)
            {
                memset(forbidden_set_matr_2[i], 0, sizeof(forbidden_set_matr_2[i]));
            }
            for (int i = 1; i <= nj; i++)
            {
                for (int j = 1; j <= nj; j++)
                {
                    memset(forbidden_set_matr_3[i][j], 0, sizeof(forbidden_set_matr_3[i][j]));
                }
            }
            vector<vector<bool>>().swap(config_matr);
            vector<vector<int>>().swap(BDP_solution);
            vector<vector<int>>().swap(solution);
            vector<vector<int>>().swap(best_solution);
            vector<int>().swap(configuration_sort_list);
        }

        freeMemoryLB4();
        FTSP_free();

    }
    return 0;
}

int sensitive_r()
{
    string out_name = ".\\sensitive-analysis-r.txt";
    vector<int> ins_dd;
    for (int i = 1; i <= 525; i++)
    {
        if (i % 75 <= 25 && i % 75 >= 1) ins_dd.push_back(i);
    }
    for (int ins_c = 0; ins_c < ins_dd.size(); ins_c++)
    {
        num_config = 3;
        cout << "instance_n = 100_" << ins_dd[ins_c] << ".alb " << endl;
        bool read_ins = instance_reader_SALBP_otto(".\\Otto dataset\\large data set_n=100\\instance_n=100_" + to_string(ins_dd[ins_c]) + ".alb");
        SALBP_trans_closure_warshall();
        SALBP_matrix_to_graph_vec();
        compute_machine_lb();
        for (int p = 1; p <= 9; p++)
        {
            string filename = ".\\FTSP dataset\\Sensitive analysis\\S_" + to_string(nj) + "_" + to_string(num_config) + "_" + to_string(p) + ".txt";
            reader_sort(filename);
            compute_C();
            configuration_order();
            heuristic_H();
            compute_forbidden_set();
            bdp();
            if (solution[0][0] > H_solution[0][0])
            {
                solution = H_solution;
            }

            if (if_generate)
            {
                best_solution = solution;
                compute_dff_lb();
                compute_LB_DW();
                LB = LB_machine > LB_1 ? LB_machine : LB_1;
                LB = LB_DW > LB ? LB_DW : LB;
                if (LB < solution[0][0] && optimalBDP != 1) ILS_solve();
                else LB = solution[0][0];
                cout << best_solution[0][0] << "\t" << LB << endl;
            }
            else
            {
                cout << "no feasible solution found!" << endl;
            }
            Writeoutlog(out_name);
            vector<vector<int>>().swap(forbidden_set_2);
            vector<vector<int>>().swap(forbidden_set_3);
            for (int i = 1; i <= nj; i++)
            {
                memset(forbidden_set_matr_2[i], 0, sizeof(forbidden_set_matr_2[i]));
            }
            for (int i = 1; i <= nj; i++)
            {
                for (int j = 1; j <= nj; j++)
                {
                    memset(forbidden_set_matr_3[i][j], 0, sizeof(forbidden_set_matr_3[i][j]));
                }
            }
            vector<vector<bool>>().swap(config_matr);
            vector<vector<int>>().swap(BDP_solution);
            vector<vector<int>>().swap(solution);
            vector<vector<int>>().swap(best_solution);
            vector<int>().swap(configuration_sort_list);
        }
        freeMemoryLB4();
        FTSP_free();
    }
    return 1;
}

int sensitive_Ck()
{
    string out_name = ".\\sensitive-analysis-Ck.txt";
    vector<int> ins_dd;
    for (int i = 1; i <= 525; i++)
    {
        if (i % 75 <= 25 && i % 75 >= 1) ins_dd.push_back(i);
    }

    for (int ins_c = 0; ins_c < ins_dd.size(); ins_c++)
    {
        num_config = 3;
        cout << "instance_n = 100_" << ins_dd[ins_c] << ".alb " << endl;
        WriteLogFile("", out_name);
        bool read_ins = instance_reader_SALBP_otto(".\\Otto dataset\\large data set_n=100\\instance_n=100_" + to_string(ins_dd[ins_c]) + ".alb");
        SALBP_trans_closure_warshall();
        SALBP_matrix_to_graph_vec();
        string filename = ".\\FTSP dataset\\S_" + to_string(nj) + "_" + to_string(num_config) + "_" + to_string(1) + ".txt";
        reader_sort(filename);
        compute_machine_lb();
        compute_C();
        configuration_order();
        compute_forbidden_set();

        int C_upper = C_k[1];
        for (int config = 2; config <= num_config; config++)
        {
            if (C_k[config] > C_upper) C_upper = C_k[config];
        }
        int C_lower = (int)ceil((double)LB_machine / (double)num_config - intTolerance);
        int C_delta = (int)floor((double)(C_upper - C_lower) / 9 + intTolerance);
        for (int step = 1; step <= 11; step++)
        {
            int C_k_value;
            if (step == 11)
            {
                C_k_value = C_lower;
            }
            else
            {
                C_k_value = C_upper - (step - 1) * C_delta;

            }
            for (int config = 1; config <= num_config; config++)
            {
                C_k[config] = C_k_value;
            }
            C_k_sum = vector<int>(num_config + 1, 0);
            for (int k = 1; k <= num_config; k++)
            {
                for (int kk = 1; kk <= k; kk++)
                {
                    C_k_sum[k] += C_k[kk];
                }
            }
            bdp();
            if (if_generate)
            {
                best_solution = solution;
                compute_dff_lb();
                compute_LB_DW();
                LB = LB_machine > LB_1 ? LB_machine : LB_1;
                LB = LB_DW > LB ? LB_DW : LB;
                if (LB < solution[0][0] && optimalBDP != 1) ILS_solve();
                else LB = solution[0][0];
                cout << best_solution[0][0] << "\t" << LB << endl;
                WriteLogFile_byEntry(to_string(step), out_name);
            }
            else
            {
                cout << "no feasible solution found!" << endl;
                break;
            }
            vector<vector<int>>().swap(BDP_solution);
            vector<vector<int>>().swap(solution);
            vector<vector<int>>().swap(best_solution);
        }
        freeMemoryLB4();
        FTSP_free();
    }
    return 1;
}

int main(void)
{
    if (version == 0) FTSP_solve();
    else if (version == 1) sensitive_r();
    else if (version == 2) sensitive_Ck();
    else exit(922);
    return 2024;
}