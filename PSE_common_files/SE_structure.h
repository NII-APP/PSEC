#ifndef SE_STRUCTURE_H
#define SE_STRUCTURE_H

#define MAXFILENAMELENGTH 15

struct SEstruct
{
	int NS; //количество граничных (внешних) узлов, переходящих на верхний уровень
	int NI; //количество исключаемых (внутренных)узлов
	int NN; //общее количество узлов
	int NSE; //количество уравнений, соответствующих NS
	int NIE; //количество уравнений, соответствующих NI
	int NNE; //общее количество узлов

	int NEL; //количество составляющих элементов
//	int neltypes; //количество типов элементов
//	int *eltypes; //номера типов элементов

	int fullnn; //количество узлов в полной модели

	int iSE; //номер суперэлемента
	int SElevel; //уровень суперэлемента

	int included_in_SE_lev; //номер уровня СЭ, в который входит данный СЭ
	int included_in_SE_num; //номер СЭ, в который входит данный СЭ


	// строка, хранящая расположение СЭ включающего данный СЭ в свой состав (требуется для расчета перемещений)
	char *PathIncludedInSE;
	char *namenumIncludedInSE;


	int NInclSEs; //количество составляющих суперэлементов
	char **PathInclSEs; //массив строк полных путей к папкам составляющих СЭ (в порядке следования в матрице индексов) (для СЭ кроме 0 уровня)
	char **namenumInclSEs;


	//флаги
	bool fl_matrix_recalc; // true - нужно выполнить повторное составление и разложение матриц
	bool fl_optrenum_performed; // true - перенумерация узлов с целью снижения объема профиля матрицы ii выполнена
	bool fl_repeated; // true - суперэлемент является повторяющимся
	bool is_matrix_calculated; // true - матрицы составлены и СЭ преобразование выполнено
	bool is_rightvect_calculated; // true - вектор нагрузок сформирован
	bool is_solvect_calculated; // true - вектор решения рассчитан
	bool is_inmem; //true - Хранение в оперативной памяти
	bool fl_PCG_internal; //true - используется сокращенное хранение STII и метод сопряженных градиентов при определении правых частей и векторов решения
	int cur_operation_type; //номер типа задания
	//2 - расчет собственных частот методом SSI

	//время
	float time_matrix;
	float time_read;
	float time_rightvect;
	float time_solvect;

	//память
	float max_required_mem;
	float min_required_mem;

	//номер задачи
	int tasknumber;

	//имена файлов
	char namenum[MAXFILENAMELENGTH];
	char name_state[MAXFILENAMELENGTH]; // файл состояния расчета для данного СЭ для использования матриц при повторном запуске
	char name_crd[MAXFILENAMELENGTH]; // матрица координат узлов - только для СЭ 0 уровня
	char name_ind[MAXFILENAMELENGTH]; // матрица индексов
	char name_mat[MAXFILENAMELENGTH]; // матрица материалов - только для СЭ 0 уровня
	char name_fix[MAXFILENAMELENGTH]; // массив закреплений - только для СЭ 0 уровня
	char name_rennodes[MAXFILENAMELENGTH]; // матрица перенумерации узлов лок-глоб
	char name_renelem[MAXFILENAMELENGTH]; // матрица перенумерации элементов лок-глоб
	char name_stii[MAXFILENAMELENGTH]; // треугольная матрица для исключаемых степеней свободы
	char name_stnz[MAXFILENAMELENGTH]; // безнулевое хранение матрицы СЭ (II или полной)
	char name_stis[MAXFILENAMELENGTH]; // прямоугольная матрица связи исключаемых и граничных степеней свободы до СЭ преобразования
	char name_stise[MAXFILENAMELENGTH]; //профиль stis в виде безнулевой схемы
	char name_stss[MAXFILENAMELENGTH]; // полностью заполненная треугольная матрица для граничных степеней свободы

	//сетевое расположение
	int cl_num;
	int cl_par_thread_num;
	char fullnetfolder[256];


	//дополнительная информация для ращличных типов расчетов

	// тип 2 - расчет собственных частот методом SSI
	double shift;

};

#endif // SE_STRUCTURE_H

/*
	char name_state[] = "state.dat"; // файл состояния расчета для данного СЭ для использования матриц при повторном запуске
	char name_crd[] = "crd.bin"; // матрица координат узлов - только для СЭ 0 уровня
	char name_ind[] = "ind.dat"; // матрица индексов
	char name_mat[] = "mat.dat"; // матрица координат узлов - только для СЭ 0 уровня
	char name_rennodes[] = "ren.dat"; // матрица перенумерации узлов лок-глоб
	char name_renelem[] = "ree.dat"; // матрица перенумерации элементов лок-глоб
	char name_stii[] = "stii.bin"; // треугольная матрица для исключаемых степеней свободы
	char name_stiie[] = "stiie.bin"; // профиль stii 
	char name_stis0[] = "stis0.bin"; // прямоугольная матрица связи исключаемых и граничных степеней свободы до СЭ преобразования
	char name_stise0[] = "stise.bin"; //профиль stis0 в виде безнулевой схемы для 0 уровня СЭ, на остальных уровнях матрица считается полностью заполненной 
	char name_stss[] = "stss.bin"; // полностью заполненная треугольная матрица для граничных степеней свободы
*/