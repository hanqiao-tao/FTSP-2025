#include "FTSP.h"

static int* apply;
static double* d_c;
static int* order;
static int remainingT;					    
static int list_pricing[NUM_JOB];		   
static double benefit_per_unit[NUM_JOB];   
static double price;					    
static double best_price;				   
static double pricing[NUM_JOB];		      
static int sol[NUM_JOB];				
static int best_sol[NUM_JOB];				
static int nj1;
static int nj2;


void initMemoryLB4()
{
    int fixed_model = kMax * (nj + 1) * (nj + 1);
    MALLOC(d_c, fixed_model, double);
    MALLOC(heads_c, fixed_model, double);
    memset(heads_c, 0, fixed_model * sizeof(double));
    MALLOC(tails_c, fixed_model, double);
    memset(tails_c, 0, fixed_model * sizeof(double));
    MALLOC(p, nj + 2, double);
    MALLOC(apply, nj + 1, int);
    MALLOC(order, nj + 1, int);

    for (int i = 1; i <= nj; i++)
    {
        p[i] = (double)(Job_Time[i]) / (double)(due_date);
    }
    nj1 = nj + 1;
    nj2 = nj1 * nj1;
}

void freeMemoryLB4()
{
    free(d_c);
    d_c = NULL;
    free(heads_c);
    heads_c = NULL;
    free(apply);
    apply = NULL;
    free(order);
    order = NULL;
    free(tails_c);
    tails_c = NULL;
    free(p);
    p = NULL;
}

inline int get_pos(int k, int eps_pos, int task)
{
    return(k * nj2 + eps_pos * nj1 + task);
}

void Sort(int list[], double* pesos, int beg, int end, int posicionInicial)
{
    int piv;
    int tmp;
    int  l, r, p;

    while (beg < end)    
    {
        l = beg; p = beg + (end - beg) / 2; r = end;
        piv = list[p];

        while (1)
        {
            while ((l <= r) && (QS_COMPARE(pesos[posicionInicial + list[l]], pesos[posicionInicial + piv]) <= 0.0)) l++;
            while ((l <= r) && (QS_COMPARE(pesos[posicionInicial + list[r]], pesos[posicionInicial + piv]) > 0.0)) r--;
            if (l > r) break;
            tmp = list[l]; list[l] = list[r]; list[r] = tmp;
            if (p == r) p = l;
            l++; r--;
        }
        list[p] = list[r]; list[r] = piv;
        r--;
        if ((r - beg) < (end - l))
        {
            Sort(list, pesos, beg, r, posicionInicial);
            beg = l;
        }
        else
        {
            Sort(list, pesos, l, end, posicionInicial);
            end = r;
        }
    }
}

void compute_data_LB4()
{
    int poisition_initial, poisition_initial_2;
    double v_cut_0, v_cut_max, f;
    for (int i = 0; i <= nj; i++)
    {
        apply[i] = 1;
    }
    for (int j = 2; j <= nj; j++)
    {
        for (int i = 1; (i < j) && (apply[j] == 1); i++)
        {
            if (Job_Time[i] == Job_Time[j])
            {
                apply[j] = 0;
            }
        }
    }
    poisition_initial = get_pos(0, 0, 0);
    d_c[poisition_initial] = 0.0;
    for (int i = 1; i <= nj; i++)
    {
        d_c[++poisition_initial] = p[i];
    }
    for (int eps_pos = 1; eps_pos <= nj; eps_pos++)
    {
        poisition_initial = get_pos(0, 0, 0);
        poisition_initial_2 = get_pos(0, eps_pos, 0);
        if (fabs(d_c[poisition_initial + eps_pos] - 0.5) < intTolerance)
        {
            v_cut_0 = 0.5;
            v_cut_max = 0.5;
        }
        else
        {
            if (d_c[poisition_initial + eps_pos] < 0.5)
            {
                v_cut_0 = d_c[poisition_initial + eps_pos];
                v_cut_max = 1.0 - d_c[poisition_initial + eps_pos];
            }
            else
            {
                v_cut_0 = 1.0 - d_c[poisition_initial + eps_pos];
                v_cut_max = d_c[poisition_initial + eps_pos];
            }
        }
        for (int i = 1; i <= nj; i++)
        {
            poisition_initial++;
            poisition_initial_2++;
            if (d_c[poisition_initial] < (v_cut_0 - intTolerance))
            {
                d_c[poisition_initial_2] = 0.0;
            }
            else
            {
                if (d_c[poisition_initial] > (v_cut_max + intTolerance))
                {
                    d_c[poisition_initial_2] = 1.0;
                }
                else
                {
                    d_c[poisition_initial_2] = d_c[poisition_initial];
                }
            }
        }
    }
    for (int k = 1; k < kMax; k++)
    {
        poisition_initial = get_pos(k, 0, 0);
        d_c[poisition_initial] = 0.0;
        for (int i = 1; i <= nj; i++)
        {
            f = p[i] * (double)(k + 1);
            if (fabs(f - round(f)) < intTolerance)
            {
                d_c[++poisition_initial] = p[i];
            }
            else
            {
                d_c[++poisition_initial] = floor(f) / (double)(k);
            }
        }
        for (int eps_pos = 1; eps_pos <= nj; eps_pos++)
        {
            poisition_initial = get_pos(k, 0, 0);
            poisition_initial_2 = get_pos(k, eps_pos, 0);
            if (fabs(d_c[poisition_initial + eps_pos] - 0.5) < intTolerance)
            {
                v_cut_0 = 0.5;
                v_cut_max = 0.5;
            }
            else
            {
                if (d_c[poisition_initial + eps_pos] < 0.5)
                {
                    v_cut_0 = d_c[poisition_initial + eps_pos];
                    v_cut_max = 1.0 - d_c[poisition_initial + eps_pos];
                }
                else
                {
                    v_cut_0 = 1.0 - d_c[poisition_initial + eps_pos];
                    v_cut_max = d_c[poisition_initial + eps_pos];
                }
            }
            for (int i = 1; i <= nj; i++)
            {
                poisition_initial++;
                poisition_initial_2++;
                if (d_c[poisition_initial] < (v_cut_0 - intTolerance))
                {
                    d_c[poisition_initial_2] = 0.0;
                }
                else
                {
                    if (d_c[poisition_initial] > (v_cut_max + intTolerance))
                    {
                        d_c[poisition_initial_2] = 1.0;
                    }
                    else
                    {
                        d_c[poisition_initial_2] = d_c[poisition_initial];
                    }
                }
            }
        }
    }
}

int compute_head(int task)
{
    int pos_order;
    if (task == 0)
    {
        pos_order = nj;
        for (int i = 1; i <= nj; i++)
        {
            order[i - 1] = nj + 1 - i;
        }
    }
    else 
    {
        pos_order = 0;
        for (int i = 1; i <= nj; i++)
        {
            if (trans_matr[i][task])
            {
                order[pos_order] = i;
                pos_order++;
            }
        }
    }
    int max_round = 0;
    for (int k_count = 0; k_count < kMax; k_count++)
    {
        for (int task_count = 0; task_count <= nj; task_count++)
        {
            if (apply[task_count] == 1)
            {
                int poisition_initial = get_pos(k_count, task_count, 0);
                if (pos_order > 1)
                {
                    Sort(order, heads_c, 0, pos_order - 1, poisition_initial);
                }
                double ac = 0.0;
                for (int i = 0; i < pos_order; i++)
                {
                    ac += d_c[poisition_initial + order[i]];
                    double tmp = ac + heads_c[poisition_initial + order[i]];
                    if (heads_c[poisition_initial + task] < tmp)
                    {
                        heads_c[poisition_initial + task] = tmp;
                    }
                }
                heads_c[poisition_initial + task] = ceil(heads_c[poisition_initial + task] - intTolerance);
                if (heads_c[poisition_initial + task] > max_round)
                {
                    max_round = heads_c[poisition_initial + task];
                }
            }
        }
    }

    for (int k_count = 0; k_count < kMax; k_count++)
    {
        for (int eps_pos = 0; eps_pos <= nj; eps_pos++)
        {
            int poisition_initial = get_pos(k_count, eps_pos, task);
            if ((heads_c[poisition_initial] < max_round) && (d_c[poisition_initial] > intTolerance))
            {
                heads_c[poisition_initial] = max_round;
            }
        }
    }
    return (0);
}

int compute_tail(int task)
{
    int pos_order;
    if (task == 0) 
    {
        pos_order = nj;
        for (int i = 1; i <= nj; i++)
        {
            order[i - 1] = i;
        }
    }
    else 
    {
        pos_order = 0;
        for (int i = 1; i <= nj; i++)
        {
            if (trans_matr[task][i])
            {
                order[pos_order] = i;
                pos_order++;
            }
        }
    }
    int max_round = 0;
    for (int k_count = 0; k_count < kMax; k_count++)
    {
        for (int task_count = 0; task_count <= nj; task_count++)
        {
            if (apply[task_count] == 1)
            {
                int poisition_initial = get_pos(k_count, task_count, 0);
                if (pos_order > 1)
                {
                    Sort(order, tails_c, 0, pos_order - 1, poisition_initial);
                }
                double ac = 0.0;
                for (int i = 0; i < pos_order; i++)
                {
                    ac += d_c[poisition_initial + order[i]];
                    double tmp = ac + tails_c[poisition_initial + order[i]];
                    if (tails_c[poisition_initial + task] < tmp)
                    {
                        tails_c[poisition_initial + task] = tmp;
                    }
                }
                tails_c[poisition_initial + task] = ceil(tails_c[poisition_initial + task] - intTolerance);
                if (tails_c[poisition_initial + task] > max_round)
                {
                    max_round = tails_c[poisition_initial + task];
                }
            }
        }
    }
    for (int k_count = 0; k_count < kMax; k_count++)
    {
        for (int eps_pos = 0; eps_pos <= nj; eps_pos++)
        {
            int poisition_initial = get_pos(k_count, eps_pos, task);
            if ((tails_c[poisition_initial] < max_round) && (d_c[poisition_initial] > intTolerance))
            {
                tails_c[poisition_initial] = max_round;
            }
        }
    }


    return (0);
}

void determine_Heads_or_Tails()
{
    clock_t start = clock();
    int pendings[NUM_JOB];
    int order_execution[NUM_JOB];

    for (int i = 1; i <= nj; i++)
    {
        order_execution[i] = nj + 1 - i;
    }
    for (int i = 1; i <= nj; i++)
    {
        pendings[i] = 0;
        for (int j = i + 1; j <= nj; j++)
        {
            if (trans_matr[i][j])
            {
                pendings[i]++;
            }
        }
    }
    int Delta = 1;
    while (Delta != 0)
    {
        Delta = 0;
        for (int i = 1; i < nj; i++)
        {
            if (pendings[order_execution[i]] > pendings[order_execution[i + 1]])
            {
                Delta = order_execution[i];
                order_execution[i] = order_execution[i + 1];
                order_execution[i + 1] = Delta;
            }
        }
    }
    for (int i = 1; i <= nj; i++)
    {
        compute_tail(order_execution[i]);
    }
    t_el += (double)(clock() - start) / CLOCKS_PER_SEC;
    compute_tail(0);

    start = clock();
    for (int i = 1; i <= nj; i++)
    {
        order_execution[i] = i;
    }
    for (int i = 1; i <= nj; i++)
    {
        pendings[i] = 0;
        for (int j = 1; j < i; j++)
        {
            if (trans_matr[j][i])
            {
                pendings[i]++;
            }
        }
    }
    Delta = 1;
    while (Delta != 0)
    {
        Delta = 0;
        for (int i = 1; i < nj; i++)
        {
            if (pendings[order_execution[i]] > pendings[order_execution[i + 1]])
            {
                Delta = order_execution[i];
                order_execution[i] = order_execution[i + 1];
                order_execution[i + 1] = Delta;
            }
        }
    }
    for (int i = 1; i <= nj; i++)
    {
        compute_head(order_execution[i]);
    }
    t_el += (double)(clock() - start) / CLOCKS_PER_SEC;
    compute_head(0);
}

bool flow(int k_value, int task_value, int LB_try)
{
    double bin_free_cap[NUM_JOB];
    double remaining;
    int task;
    if (LB_try > nj)
    {
        return true;
    }
    int position_initial = get_pos(k_value, task_value, 0);
    for (int i = 1; i <= LB_try; i++)
    {
        bin_free_cap[i] = 1.0;
    }
    for (int i = 1; i <= nj; i++)
    {
        task = order[i];
        remaining = d_c[position_initial + task];
        for (int j = E[task]; (j <= L[task]) && (remaining > intTolerance); j++)
        {
            if (remaining < bin_free_cap[j])
            {
                bin_free_cap[j] -= remaining;
                remaining = 0.0;
            }
            else
            {
                remaining -= bin_free_cap[j];
                bin_free_cap[j] = 0.0;
            }
        }
        if (remaining > intTolerance)
        {
            return false;
        }
    }
    return true;
}

int determine_flow(int LB)
{
    int LB_try = LB;
    for (int i = 1; i <= nj; i++)
    {
        order[i] = i;
    }
    for (int i = 1; i <= nj; i++)
    {
        E[i] = 1 + (int)(heads_c[i]);
    }
repetir:
    for (int i = 1; i <= nj; i++)
    {
        L[i] = LB_try - (int)(tails_c[i]);
    }
    int Delta = 1;
    while (Delta > 0)
    {
        Delta = 0;
        for (int i = 1; i < nj; i++)
        {
            if (L[order[i]] > L[order[i + 1]])
            {
                Delta = order[i];
                order[i] = order[i + 1];
                order[i + 1] = Delta;
            }
            else
            {
                if ((L[order[i]] == L[order[i + 1]]) && (E[order[i]] > E[order[i + 1]]))
                {
                    Delta = order[i];
                    order[i] = order[i + 1];
                    order[i + 1] = Delta;
                }
            }
        }
    }
    for (int k_count = 0; k_count < kMax; k_count++)
    {
        for (int task_count = 0; task_count <= nj; task_count++)
        {
            if (flow(k_count, task_count, LB_try) == false)
            {
                LB_try++;
                goto repetir;
            }
        }
    }
    return LB_try;
}

void compute_machine_lb() 
{
    clock_t start = clock();
    initMemoryLB4();
    compute_data_LB4();
    t_el += (double)(clock() - start) / CLOCKS_PER_SEC;
    determine_Heads_or_Tails();
    int LB_0 = (int)(tails_c[0]);
    if (LB_0 < (int)(heads_c[0]))
    {
        LB_0 = (int)(heads_c[0]);
    }
    LB_machine = determine_flow(LB_0);
    t_machine = (double)(clock() - start) / CLOCKS_PER_SEC;
}

void compute_dff_lbmodel(int kk, int eps_pos)
{
    IloEnv env;
    try {
        IloModel model(env);

        NumVarMatrix x(env, nj + 1);
        for (int i = 0; i <= nj; i++)
        {
            x[i] = IloNumVarArray(env);
            for (int j = 0; j <= solution[0][0]; j++)
            {
                string name = "x" + to_string(i) + to_string(j);
                x[i].add(IloNumVar(env, 0, 1, ILOFLOAT, name.c_str()));
            }
            x[i][0].setBounds(0, 0);
        }
        for (int j = 0; j <= solution[0][0]; j++)
        {
            x[0][j].setBounds(0, 0);
        }

        NumVarMatrix y(env, solution[0][0] + 1);
        for (int j = 0; j <= solution[0][0]; j++)
        {
            y[j] = IloNumVarArray(env);
            for (int k = 0; k <= num_config; k++)
            {
                string name = "y" + to_string(j) + to_string(k);
                y[j].add(IloNumVar(env, 0, 1, ILOFLOAT, name.c_str()));
            }
            y[j][0].setBounds(0, 0);
        }
        for (int k = 0; k <= num_config; k++)
        {
            y[0][k].setBounds(0, 0);
        }

        IloExpr obj(env);
        for (int j = 0; j <= solution[0][0]; j++)
        {
            obj += IloSum(y[j]);
        }
        model.add(IloMinimize(env, obj));
        obj.end();

        IloRangeArray constraints(env);
        for (int i = 1; i <= nj; i++)
        {
            constraints.add(IloSum(x[i]) == 1);
        }
        for (int j = 1; j <= solution[0][0]; j++)
        {
            constraints.add(IloSum(y[j]) <= 1);
        }

        for (int j = 1; j <= solution[0][0]; j++)
        {
            IloExpr sum_x_j(env);
            for (int i = 1; i <= nj; i++)
            {
                sum_x_j += (d_c[get_pos(kk, eps_pos, i)] * x[i][j]);
            }
            sum_x_j -= IloSum(y[j]);
            constraints.add(sum_x_j <= 0);
            sum_x_j.end();
        }

        for (int pre = 0; pre < pre_rela.size(); pre++)
        {
            int j1 = pre_rela[pre].first;
            int j2 = pre_rela[pre].second;
            IloExpr prec(env);
            for (int i = 1; i <= solution[0][0]; i++)
            {
                prec += i * x[j2][i];
                prec -= i * x[j1][i];
            }
            constraints.add(prec >= 1);
            prec.end();
        }

        for (int i = 1; i <= nj; i++)
        {
            for (int j = 1; j <= solution[0][0]; j++)
            {
                IloExpr rhs(env);
                for (int k = 1; k <= num_config; k++)
                {
                    if (config_matr[k][i] == 1)
                    {
                        rhs += y[j][k];
                    }
                }
                constraints.add(x[i][j] - rhs <= 0);
                rhs.end();
            }
        }

        for (int k = 1; k <= num_config; k++)
        {
            IloExpr lhs(env);
            for (int j = 1; j <= solution[0][0]; j++)
            {
                lhs += y[j][k];
            }
            constraints.add(lhs - C_k[k] <= 0);
            lhs.end();
        }

        for (int i = 1; i <= nj; i++)
        {
            for (int j = 1; j < E[i]; j++)
            {
                x[i][j].setBounds(0, 0);
            }
            for (int j = solution[0][0]; j > L[i]; j--)
            {
                x[i][j].setBounds(0, 0);
            }
        }
        model.add(constraints);
        constraints.end();

        IloCplex cplex(model);
        cplex.setParam(IloCplex::Param::Threads, 1);
        cplex.setOut(env.getNullStream());
        cplex.solve();
        if (cplex.getStatus() == IloAlgorithm::Optimal)
            LB_1 = LB_1 > (int)ceil(cplex.getObjValue() - intTolerance) ? LB_1 : (int)ceil(cplex.getObjValue() - intTolerance);
        else LB_1 = 0;
    }
    catch (IloException& e)
    {
        cerr << "Concert exception caught: " << e << endl;
    }
    catch (...)
    {
        cerr << "Unknown exception caught" << endl;
    }
    env.end();
}

void compute_dff_lb()
{
    clock_t start = clock();
    LB_1 = 0;
    if (nj > 100) compute_dff_lbmodel(0, 0);
    else
    {
        for (int k_count = 0; k_count < kMax; k_count++)
        {
            compute_dff_lbmodel(k_count, 0);
        }
    }
    t_lb_1 = (double)(clock() - start) / CLOCKS_PER_SEC;
}

double bound_pricing(int rest_time, int position)
{
    int i, task;
    double ac = 0.0;
    for (i = position; i <= nj; i++)
    {
        task = list_pricing[i];
        if (rest_time > Job_Time[task])
        {
            rest_time -= Job_Time[task];
            ac += pricing[task];
        }
        else
        {
            ac += pricing[task] * (double)(rest_time) / (double)(Job_Time[task]);
            return ac;
        }
    }
    return ac;
}

int knapsack_branch(int position, int num_select)
{
    int task = list_pricing[position];
    int i, j;
    if (Job_Time[task] > remainingT) return 0;
    for (i = 0; i < num_select; i++)
    {
        if (trans_matr[task][sol[i]]) return 0;
        if (trans_matr[sol[i]][task]) return 0;
    }
    for (i = 0; i < num_select; i++)
    {
        if (forbidden_set_matr_2[task][sol[i]]) return 0;
    }
    for (i = 0; i < num_select; i++)
    {
        for (j = i + 1; j < num_select; j++)
        {
            if (forbidden_set_matr_3[sol[i]][sol[j]][task]) return 0;
        }
    }
    int eMin = E[task];
    int lMax = L[task];
    for (i = 0; i < num_select; i++) if (eMin < E[sol[i]]) eMin = E[sol[i]];
    for (i = 0; i < num_select; i++) if (lMax > L[sol[i]]) lMax = L[sol[i]];
    if (eMin > lMax) return(0);
    if (price + pricing[task] + bound_pricing(remainingT - Job_Time[task], position + 1) < best_price) return 0;
    sol[num_select] = task;
    remainingT -= Job_Time[task];
    price += pricing[task];
    num_select++;
    if (price > best_price)
    {
        memset(best_sol, 0, sizeof(best_sol));
        for (int i = 0; i < num_select; i++)
        {
            best_sol[sol[i]] = 1;
        }
        best_price = price;
    }
    for (i = position + 1; i <= nj; i++)
    {
        if (pricing[list_pricing[i]] > intTolerance)
        {
            knapsack_branch(i, num_select);
        }
    }

    num_select--;
    price -= pricing[task];
    remainingT += Job_Time[task];
    sol[num_select] = 0;

    return 0;

}

void compute_forbidden_set()
{
    forbidden_set_2.reserve(nj * (nj - 1) / 2);
    vector<vector<int>> nonforbidden_set_2;
    nonforbidden_set_2.reserve(nj * (nj - 1) / 2);
    for (int j1 = 1; j1 <= nj; j1++)
    {
        for (int j2 = j1 + 1; j2 <= nj; j2++)
        {
            bool flag = true;
            for (int k = 1; flag && k <= num_config; k++)
            {
                if (config_matr[k][j1] && config_matr[k][j2])
                {
                    flag = false;
                }
            }
            vector<int> tmp = { j1,j2 };
            if (flag)
            {
                forbidden_set_2.push_back(tmp);
                forbidden_set_matr_2[j1][j2] = 1;
                forbidden_set_matr_2[j2][j1] = 1;
            }
            else
            {
                nonforbidden_set_2.push_back(tmp);
            }
        }
    }

    if (num_config > 2)
    {
        forbidden_set_3.reserve(nj * (nj - 1) * (nj - 2) / 6);
        for (int i = 0; i < nonforbidden_set_2.size(); i++)
        {
            int j1 = nonforbidden_set_2[i][0];
            int j2 = nonforbidden_set_2[i][1];
            for (int j = j2 + 1; j <= nj; j++)
            {
                bool flag = true;
                for (int k = 1; flag && k <= num_config; k++)
                {
                    if ((config_matr[k][j1] && config_matr[k][j2] && config_matr[k][j]) || forbidden_set_matr_2[j1][j] == 1 || forbidden_set_matr_2[j2][j] == 1)
                    {
                        flag = false;
                    }
                }
                if (flag)
                {
                    vector<int>tmp = { j1,j2,j };
                    forbidden_set_3.push_back(tmp);
                }
            }
        }
    }

    for (int i = 0; i < forbidden_set_3.size(); i++)
    {
        int j1 = forbidden_set_3[i][0];
        int j2 = forbidden_set_3[i][1];
        int j3 = forbidden_set_3[i][2];
        forbidden_set_matr_3[j1][j2][j3] = 1;
        forbidden_set_matr_3[j1][j3][j2] = 1;
        forbidden_set_matr_3[j2][j1][j3] = 1;
        forbidden_set_matr_3[j2][j3][j1] = 1;
        forbidden_set_matr_3[j3][j1][j2] = 1;
        forbidden_set_matr_3[j3][j2][j1] = 1;
    }
}

void pricing_solve()
{
    clock_t start = clock();
    int i, num_select;
    int change = 1;
    price = 0.0;
    best_price = 1.0;
    for (i = 1; i <= nj; i++)
    {
        benefit_per_unit[i] = pricing[i] / (double)Job_Time[i];
        list_pricing[i] = i;
    }
    while (change >= 0)
    {
        change = -1;
        for (int i = 1; i < nj; i++)
        {
            if (benefit_per_unit[list_pricing[i]] < benefit_per_unit[list_pricing[i + 1]])
            {
                change = list_pricing[i];
                list_pricing[i] = list_pricing[i + 1];
                list_pricing[i + 1] = change;
            }
        }
    }
    for (i = 1; i <= nj; i++)
    {
        sol[i] = 0;
        best_sol[i] = 0;
    }
    remainingT = due_date;
    num_select = 0;
    for (i = 1; i <= nj; i++)
    {
        if (pricing[list_pricing[i]] > intTolerance)
        {
            knapsack_branch(i, num_select);
        }
    }
    t_pricing += (double)(clock() - start) / CLOCKS_PER_SEC;
}

void compute_LB_DW()
{
    clock_t start = clock();
    IloEnv env;
    try
    {
        IloModel model(env);
        IloObjective obj = IloMinimize(env);
        model.add(obj);
        IloNumVarArray z_q(env);
        IloRangeArray col_constraints(env);
        IloNumArray pi(env, best_solution[0][0] + 1); 
        IloExpr empty(env);
        col_constraints.add(empty >= 0);
        empty -= 1;
        for (int i = 1; i <= nj; i++)
        {
            col_constraints.add(empty == 0);
        }
        empty.end();
        model.add(col_constraints);
        IloNumArray new_column(env, nj + 1);
        new_column[0] = 0;
        for (int j = 1; j <= best_solution[0][0]; j++)
        {
            for (int i = 1; i <= nj; i++)
            {
                if (best_solution[0][i] == j)
                {
                    new_column[i] = 1;
                }
                else
                {
                    new_column[i] = 0;
                }
            }
            z_q.add(IloNumVar(obj(1) + col_constraints(new_column), 0, 1, ILOFLOAT));
        }
        new_column.end();
        IloCplex cplex(model);
        cplex.setParam(IloCplex::Param::Threads, 1);
        cplex.setOut(env.getNullStream());
        cplex.solve();

        cplex.getDuals(pi, col_constraints);
        bool col_gen = false;
        pricing[0] = 0.0;
        for (int i = 1; i <= nj; i++)
        {
            pricing[i] = pi[i];
        }
        pricing_solve();
        if (best_price - 1 > intTolerance)
        {
            IloNumArray new_column(env, nj + 1);
            new_column[0] = 0;
            for (int i = 1; i <= nj; i++)
            {
                if (best_sol[i] == 1)
                {
                    new_column[i] = 1;
                }
                else new_column[i] = 0;
            }
            z_q.add(IloNumVar(obj(1) + col_constraints(new_column), 0, 1, ILOFLOAT));
            new_column.end();
            col_gen = true;
        }

        while (col_gen)
        {
            col_gen = false;
            cplex.solve();
            cplex.getDuals(pi, col_constraints);
            for (int i = 1; i <= nj; i++)
            {
                pricing[i] = pi[i];
            }
            pricing_solve();
            if (best_price - 1 > intTolerance)
            {
                IloNumArray new_column(env, nj + 1);
                new_column[0] = 0;
                for (int i = 1; i <= nj; i++)
                {
                    if (best_sol[i] == 1)
                    {
                        new_column[i] = 1;
                    }
                    else new_column[i] = 0;
                }
                z_q.add(IloNumVar(obj(1) + col_constraints(new_column), 0, 1, ILOFLOAT));
                new_column.end();
                col_gen = true;
            }
        }
        z_q.end();
        LB_DW = (int)(ceil(cplex.getObjValue() - intTolerance));
        LB = LB_DW > LB ? LB_DW : LB;
    }
    catch (IloException& e)
    {
        cerr << "Concert exception caught: " << e << endl;
    }
    catch (...)
    {
        cerr << "Unknown exception caught" << endl;
    }
    env.end();
    t_dw = t_el + (double)(clock() - start) / CLOCKS_PER_SEC;
}