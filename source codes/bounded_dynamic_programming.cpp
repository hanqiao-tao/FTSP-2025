#include "FTSP.h"

static char* degrees;
static int* eligible;
static int* tasks;
static int* list_do;                    
static double* tails;                  
static int max_loads;
static int n_full_loads;
static int windowsWidth = 1000;        
static int maxTransit = 1000;
static int maxTries = 500;
static int station;
static int stationAssignment[NUM_JOB];
static int workingDegrees[NUM_JOB];
static int stationConfiguration[NUM_JOB];

static Data elements;
static Data State_enumera[600000];
static priority_queue<Data, std::vector<Data>, compare_pq > Heap_Old;
static priority_queue<Data, std::vector<Data>, compare_pq > Heap_Enumera;
static priority_queue<Data, std::vector<Data>, compare_pq > Heap_New;
static set<Data, compare_set > State;

static int* Bipartite[NUM_JOB];         
static int* match;					   
static bool* vis;					    


bool dfs(int u)
{
    int i;
    for (i = 1; i <= Bipartite[u][0]; i++) 
    {
        int v = Bipartite[u][i];
        if (!vis[v]) {
            vis[v] = true;
            if (match[v] == -1 || dfs(match[v])) 
            {
                match[v] = u;
                return true;
            }
        }
    }
    return false;
}

bool hungarian()
{
    int res = 0;

    for (int i = 0; i <= nj; i++)
    {
        MALLOC(Bipartite[i], C_k_sum[num_config] + 1, int);
        Bipartite[i][0] = 0;
    }
    vector<vector<int>> feasible_config(station + 1, vector<int>(num_config + 1, 1));
    for (int i = 1; i <= nj; i++)
    {
        for (int k = 1; k <= num_config; k++)
        {
            if (!config_matr[k][i])
            {
                feasible_config[stationAssignment[i]][k] = 0;
            }
        }
    }
    for (int u = 1; u <= station; u++)
    {
        for (int k = 1; k <= num_config; k++)
        {
            if (feasible_config[u][k] == 1)
            {
                for (int r = 1; r <= C_k[k]; r++)
                {
                    Bipartite[u][0]++;
                    Bipartite[u][Bipartite[u][0]] = C_k_sum[k - 1] + r;
                }
            }
        }
    }

    MALLOC(match, C_k_sum[num_config] + 1, int);
    MALLOC(vis, C_k_sum[num_config] + 1, bool);
    memset(match, -1, (C_k_sum[num_config] + 1) * sizeof(int));
    for (int u = 1; u <= station; u++)
    {
        memset(vis, false, (C_k_sum[num_config] + 1) * sizeof(bool));
        if (dfs(u)) res++;
    }
    for (int i = 0; i <= nj; i++)
    {
        free(Bipartite[i]);
        Bipartite[i] = NULL;
    }

    free(vis);
    vis = NULL;
    if (res == station) return true;
    return false;
}

void init_bounds()
{
    MALLOC(list_do, nj + 1, int);
    MALLOC(tails, nj + 1, double);
    for (int i = 1; i <= nj; i++) list_do[i] = i;
    int delta = 1;
    while (delta > 0)
    {
        delta = 0;
        for (int i = 1; i <= nj - 1; i++)
        {
            if (tails_c[list_do[i]] < tails_c[list_do[i + 1]])
            {
                delta = list_do[i];
                list_do[i] = list_do[i + 1];
                list_do[i + 1] = delta;
            }
        }
    }
}

void free_bounds()
{
    free(list_do);
    list_do = NULL;
    free(tails);
    tails = NULL;
}

int bound4(int ne, char* deg)
{
    double ac = 0.0;
    double bestbound = 0.0;
    int bestSP = 0;
    int n_pendings = 0;
    int t_pendings = 0;
    int sum_successor_pendings = 0;

    for (int i = 1; i <= nj; i++)
    {
        if (deg[list_do[i]] >= 0)
        {
            if (bestSP < shortest_path[list_do[i]]) bestSP = shortest_path[list_do[i]];
            t_pendings += Job_Time[list_do[i]];
            sum_successor_pendings += all_sucs[list_do[i]].size();
            n_pendings++;
            ac += p[list_do[i]];
            if ((ac + tails_c[list_do[i]]) > bestbound)
            {
                bestbound = ac + tails_c[list_do[i]];
            }
        }
    }
    elements.weight = (bestbound + ne) * 1000000.0 + bestSP * 1000.0 - n_pendings;
    return((int)(ceil(bestbound - intTolerance)) + ne);
}

void initialize_bdp()
{
    MALLOC(degrees, nj + 1, char);
    MALLOC(eligible, nj + 1, int);
    MALLOC(tasks, nj + 1, int);
}

bool saveBestSolution()
{
    int i, j, estacion, esOK;
    memset(workingDegrees, -1, sizeof(workingDegrees));
    memset(stationAssignment, -1, sizeof(stationAssignment));
    memset(stationConfiguration, -1, sizeof(stationConfiguration));

    for (i = station; i > 0; i--)
    {
        for (j = 1; j <= nj; j++)
        {
            if (elements.degrees[j] == -1 && elements.previous_ptr->degrees[j] == 0)
            {
                stationAssignment[j] = i;
            }
        }
        elements = *elements.previous_ptr;
    }
    bool flag = hungarian();
    if (flag)
    {
        for (int i = 1; i <= C_k_sum[num_config]; i++)
        {
            if (match[i] != -1)
            {
                for (int k = 0; k < num_config; k++)
                {
                    if (i > C_k_sum[k] && i <= C_k_sum[k + 1])
                    {
                        stationConfiguration[match[i]] = k + 1;
                        break;
                    }
                }
            }
        }
    }
    free(match);
    match = NULL;
    return flag;
}

void test_solution_complete()
{
    int t_unassigned = 0;
    int i;
    for (i = 1; (i <= nj) && (t_unassigned == 0); i++)
    {
        if (degrees[i] >= 0) t_unassigned += Job_Time[i];
    }
    memcpy(elements.degrees, degrees, nj + 1);
    if ((t_unassigned == 0) && (station < UB_BDP))
    {

        bool flag = saveBestSolution();
        if (flag)
        {
            UB_BDP = station;
            BDP_solution.push_back(vector<int>(nj + 1));
            BDP_solution[0][0] = UB_BDP;
            BDP_solution.push_back(vector<int>(UB_BDP + 1));
            {
                for (i = 1; i <= nj; i++)
                {
                    BDP_solution[0][i] = stationAssignment[i];
                }
                for (i = 1; i <= UB_BDP; i++)
                {
                    BDP_solution[1][i] = stationConfiguration[i];
                }
            }
        }
    }
}

int task_solitary(int n_eligible)
{
    int result, depth, i, j, jj;
    for (i = 1; i <= n_eligible; i++)
    {
        int task = eligible[i];
        result = 1; 
        for (j = 1; (j <= nj) && (result == 1); j++)
        {
            if (degrees[j] >= 0 && (task != j) && (trans_matr[task][j] == false) && (trans_matr[j][task] == false) && Job_Time[task] + Job_Time[j] <= due_date && !forbidden_set_matr_2[task][j])
            {
                result = 0; 
            }
        }
        if (result == 1)
        {
            degrees[task] = -1;
            int n_sub_eligible = n_eligible;
            for (jj = 0; jj < direct_sucs[task].size(); jj++)
            {
                j = direct_sucs[task][jj];
                degrees[j]--;
            }
            test_solution_complete();
            if ((bound4(station, degrees)) >= UB_BDP) return(1);
            Heap_Enumera.push(elements);
            return(1);
        }
    }
    return (0);
}

void gen_loads_bdp(int depth, int remaining_time, int start, int n_eligible, vector<int> config_fea)
{
    int full_load = 1;
    int n_sub_eligible, i, ii, j, jj;

    for (ii = start; ii <= n_eligible; ii++)
    {
        i = eligible[ii];
        if (Job_Time[i] <= remaining_time)
        {
            vector<int> config_fea_old = config_fea;
            int if_config = 0;
            for (int k = 1; k <= num_config; k++)
            {
                if (config_matr[k][i] && config_fea[k])
                {
                    if_config++;
                }
                else
                {
                    config_fea[k] = 0;
                }
            }
            if (if_config > 0)
            {
                tasks[depth] = i;
                degrees[i] = -1;
                n_sub_eligible = n_eligible;
                full_load = 0;
                for (jj = 0; jj < direct_sucs[i].size(); jj++)
                {
                    j = direct_sucs[i][jj];
                    degrees[j]--;
                }
                gen_loads_bdp(depth + 1, remaining_time - Job_Time[i], ii + 1, n_sub_eligible, config_fea);
                if (nj >= 100 && (n_full_loads >= max_loads))
                {
                    optimalBDP = 0;
                    return;
                }
                tasks[depth] = -1;
                degrees[i] = 0;
                config_fea = config_fea_old;
                for (jj = 0; jj < direct_sucs[i].size(); jj++)
                {
                    j = direct_sucs[i][jj];
                    degrees[j]++;
                }
            }
        }
    }
    for (ii = 1; (ii <= n_eligible) && (full_load == 1); ii++)
    {
        i = eligible[ii];
        if ((degrees[i] == 0) && (remaining_time >= Job_Time[i]))
        {
            for (int k = 1; (k <= num_config) && (full_load == 1); k++)
            {
                if (config_matr[k][i] && config_fea[k])
                {
                    full_load = 0;
                }
            }
        }
    }

    if (full_load == 1)
    {
        n_full_loads++;
        if (depth == 1) return;
        test_solution_complete();
        ii = 1;
        if ((bound4(station, degrees)) >= UB_BDP) ii = 0;
        if (ii == 1)
        {
            if (Heap_Enumera.size() < maxTransit)
            {
                Heap_Enumera.push(elements);
            }
            else
            {
                optimalBDP = 0;
                if (elements.weight < Heap_Enumera.top().weight)
                {
                    Heap_Enumera.pop();
                    Heap_Enumera.push(elements);
                }
            }
        }
    }
}

void free_bdp()
{
    free(degrees);
    degrees = NULL;
    free(eligible);
    eligible = NULL;
    free(tasks);
    tasks = NULL;
    while (!Heap_New.empty()) Heap_New.pop();
    while (!Heap_Old.empty()) Heap_Old.pop();
    while (!Heap_Enumera.empty()) Heap_Enumera.pop();
    State.clear();
    set<Data, compare_set >().swap(State);
    memset(State_enumera, 0, sizeof(State_enumera));
}

void bdp()
{
    init_bounds();
    initialize_bdp();
    optimalBDP = 1;
    int i, j, count, index, step, n_eligible, pos_count_1;
    State.clear();
    max_loads = maxTries;
    clock_t  start_time = clock();
    for (i = 0; i <= nj; i++) elements.degrees[i] = 0x00;
    UB_BDP = nj;
    for (i = 1; i <= nj; i++)
    {
        count = 0;
        for (j = 1; j < i; j++)
        {
            if (adj_matr[j][i] == 1) count++;
        }
        degrees[i] = count;
    }
    memcpy(elements.degrees, degrees, nj + 1);
    pos_count_1 = 0;
    Heap_New.push(elements);
    for (step = 1; step <= UB_BDP; step++)
    {
        station = step;
        if (Heap_New.empty()) break;
        while (!Heap_New.empty())
        {
            elements = Heap_New.top();
            State_enumera[pos_count_1] = elements;
            elements.previous_ptr = &State_enumera[pos_count_1];
            pos_count_1++;
            Heap_Old.push(elements);
            Heap_New.pop();
        }
        while (!Heap_Old.empty())
        {
            elements = Heap_Old.top();
            Heap_Old.pop();
            for (i = 1; i <= nj; i++) degrees[i] = elements.degrees[i];
            n_eligible = 0;
            n_full_loads = 0;
            for (i = 1; i <= nj; i++)
            {
                if (degrees[i] == 0)
                {
                    eligible[++n_eligible] = i;
                }
            }

            if (task_solitary(n_eligible) == 0)
            {
                vector<int>config_fea(num_config + 1, 1);
                gen_loads_bdp(1, due_date, 1, n_eligible, config_fea);
            }
            while (!Heap_Enumera.empty())
            {
                elements = Heap_Enumera.top();
                Heap_Enumera.pop();
                if (State.find(elements) == State.end())
                {
                    if (Heap_New.size() < windowsWidth)
                    {
                        State.insert(elements);
                        Heap_New.push(elements);
                    }
                    else
                    {
                        optimalBDP = 0;
                        if (Heap_New.top().weight > elements.weight)
                        {
                            State.insert(elements);
                            Heap_New.pop();
                            Heap_New.push(elements);
                        }
                    }
                }
            }
        }
        State.clear();
    }
    t_BDP = (double)(clock() - start_time) / CLOCKS_PER_SEC;

    for (int i = 1; i <= nj; i++)
    {
        E[i] = 1 + floor(heads_c[i] + intTolerance);
        L[i] = UB_BDP - floor(tails_c[i] + intTolerance);
    }
    free_bdp();
    free_bounds();
    if (BDP_solution.size() > 0)
    {
        bool flag = true;
        for (int i = 1; flag && i <= nj; i++)
        {
            if (BDP_solution[0][i] == -1)
            {
                if_generate = 0;
                flag = false;
            }
        }
        if (flag) solution = BDP_solution;
    }
    else if_generate = 0;
}