#ifndef SE_STRUCTURE_H
#define SE_STRUCTURE_H

#define MAXFILENAMELENGTH 15

struct SEstruct
{
	int NS; //���������� ��������� (�������) �����, ����������� �� ������� �������
	int NI; //���������� ����������� (����������)�����
	int NN; //����� ���������� �����
	int NSE; //���������� ���������, ��������������� NS
	int NIE; //���������� ���������, ��������������� NI
	int NNE; //����� ���������� �����

	int NEL; //���������� ������������ ���������
//	int neltypes; //���������� ����� ���������
//	int *eltypes; //������ ����� ���������

	int fullnn; //���������� ����� � ������ ������

	int iSE; //����� �������������
	int SElevel; //������� �������������

	int included_in_SE_lev; //����� ������ ��, � ������� ������ ������ ��
	int included_in_SE_num; //����� ��, � ������� ������ ������ ��


	// ������, �������� ������������ �� ����������� ������ �� � ���� ������ (��������� ��� ������� �����������)
	char *PathIncludedInSE;
	char *namenumIncludedInSE;


	int NInclSEs; //���������� ������������ ��������������
	char **PathInclSEs; //������ ����� ������ ����� � ������ ������������ �� (� ������� ���������� � ������� ��������) (��� �� ����� 0 ������)
	char **namenumInclSEs;


	//�����
	bool fl_matrix_recalc; // true - ����� ��������� ��������� ����������� � ���������� ������
	bool fl_optrenum_performed; // true - ������������� ����� � ����� �������� ������ ������� ������� ii ���������
	bool fl_repeated; // true - ������������ �������� �������������
	bool is_matrix_calculated; // true - ������� ���������� � �� �������������� ���������
	bool is_rightvect_calculated; // true - ������ �������� �����������
	bool is_solvect_calculated; // true - ������ ������� ���������
	bool is_inmem; //true - �������� � ����������� ������
	bool fl_PCG_internal; //true - ������������ ����������� �������� STII � ����� ����������� ���������� ��� ����������� ������ ������ � �������� �������
	int cur_operation_type; //����� ���� �������
	//2 - ������ ����������� ������ ������� SSI

	//�����
	float time_matrix;
	float time_read;
	float time_rightvect;
	float time_solvect;

	//������
	float max_required_mem;
	float min_required_mem;

	//����� ������
	int tasknumber;

	//����� ������
	char namenum[MAXFILENAMELENGTH];
	char name_state[MAXFILENAMELENGTH]; // ���� ��������� ������� ��� ������� �� ��� ������������� ������ ��� ��������� �������
	char name_crd[MAXFILENAMELENGTH]; // ������� ��������� ����� - ������ ��� �� 0 ������
	char name_ind[MAXFILENAMELENGTH]; // ������� ��������
	char name_mat[MAXFILENAMELENGTH]; // ������� ���������� - ������ ��� �� 0 ������
	char name_fix[MAXFILENAMELENGTH]; // ������ ����������� - ������ ��� �� 0 ������
	char name_rennodes[MAXFILENAMELENGTH]; // ������� ������������� ����� ���-����
	char name_renelem[MAXFILENAMELENGTH]; // ������� ������������� ��������� ���-����
	char name_stii[MAXFILENAMELENGTH]; // ����������� ������� ��� ����������� �������� �������
	char name_stnz[MAXFILENAMELENGTH]; // ���������� �������� ������� �� (II ��� ������)
	char name_stis[MAXFILENAMELENGTH]; // ������������� ������� ����� ����������� � ��������� �������� ������� �� �� ��������������
	char name_stise[MAXFILENAMELENGTH]; //������� stis � ���� ���������� �����
	char name_stss[MAXFILENAMELENGTH]; // ��������� ����������� ����������� ������� ��� ��������� �������� �������

	//������� ������������
	int cl_num;
	int cl_par_thread_num;
	char fullnetfolder[256];


	//�������������� ���������� ��� ��������� ����� ��������

	// ��� 2 - ������ ����������� ������ ������� SSI
	double shift;

};

#endif // SE_STRUCTURE_H

/*
	char name_state[] = "state.dat"; // ���� ��������� ������� ��� ������� �� ��� ������������� ������ ��� ��������� �������
	char name_crd[] = "crd.bin"; // ������� ��������� ����� - ������ ��� �� 0 ������
	char name_ind[] = "ind.dat"; // ������� ��������
	char name_mat[] = "mat.dat"; // ������� ��������� ����� - ������ ��� �� 0 ������
	char name_rennodes[] = "ren.dat"; // ������� ������������� ����� ���-����
	char name_renelem[] = "ree.dat"; // ������� ������������� ��������� ���-����
	char name_stii[] = "stii.bin"; // ����������� ������� ��� ����������� �������� �������
	char name_stiie[] = "stiie.bin"; // ������� stii 
	char name_stis0[] = "stis0.bin"; // ������������� ������� ����� ����������� � ��������� �������� ������� �� �� ��������������
	char name_stise0[] = "stise.bin"; //������� stis0 � ���� ���������� ����� ��� 0 ������ ��, �� ��������� ������� ������� ��������� ��������� ����������� 
	char name_stss[] = "stss.bin"; // ��������� ����������� ����������� ������� ��� ��������� �������� �������
*/