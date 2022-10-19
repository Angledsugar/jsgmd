"""
Author Kcrong
"""
from random import uniform, randint

from matplotlib.pyplot import plot, show, xlim, ylim, xlabel, ylabel
from numpy import mean

# According to PEP8, do not assign lambda
def rand(x, y): return int(uniform(x, y))

# 유전자 갯수
_gen_len = 30

# 유전자 내 최소 및 최댓값
_gen_min_value = 1
_gen_max_value = 30

# 우월 유전자 보존 갯수
GOOD_DNA_CNT = 6

# 돌연변이 확률은 fitness 와 반비례 한다.
# fitness 가 높을 수록, 돌연변이 확률이 적어진다.
MUTATION_PROBABILITY = 10


class Generation:
    cnt = 0

    def __init__(self, dna_list):
        Generation.cnt += 1
        self.generation_level = Generation.cnt
        self.DNA_list = dna_list
        self.select_list = self.make_select_list()

    def __repr__(self):
        return "<Generation level %d>" % self.generation_level

    def make_select_list(self):
        """
        dna fitness 만큼의 갯수를 가지는 dna 리스트
        dna1.fitness = 2,
        dna2.fitness = 3, then
        return [dna1, dna1, dna2, dna2, dna2]
        """

        tmp_list = list()

        for dna in self.DNA_list:
            tmp_list += [dna for _ in range(dna.fitness)]

        return tmp_list

    def make_child(self):
        """
        :return: Child Gene Object
        """

        # 돌연변이(mutation)
        # fitness가 높아질수록, 돌연변이 확률이 적어짐
        if rand(0, self.fitness * MUTATION_PROBABILITY) == 0:
            gen_mutate = []
            ran_num = randint(_gen_min_value, _gen_max_value)

            for i in range(_gen_len):
                while ran_num in gen_mutate:
                    ran_num = randint(_gen_min_value, _gen_max_value)
                gen_mutate.append(ran_num)

            return DNA(gen_mutate)

        # 부모를 select_list 를 이용해 정함.
        # 부모로 선출될 확률은 fitness 과 비례한다. 
        # 룰렛 휠 선택 방식 : 적합도가 높을수록 부모 해가 선택될 확률이 높아지는 방식
        parents = tuple(self.select_list[rand(0, len(self.select_list))] for _ in range(2))

        # 자식 유전자
        gene_data = list()

        # 유전자 정보 길이
        gene_len = len(parents[0].gene_data)

        # 각 교차 포인트를 정한다.
        # rand 함수의 반환이 float 형이므로, 소수점을 버리기 위해 int() 형변한을 해준다.
        switch_point = (rand(1, gene_len // 2), rand(gene_len // 2, gene_len))

        # 처음 자식이 받는 유전자는 parent1
        # 다만 교차 포인트에 다다르면, 다른 parent 유전자 정보를 받아오기 시작한다. (parent = parent2)
        parent = parents[0]

        for _ in range(gene_len):
            # 자식 유전자 정보는 부모 유전자에서 받아온다
            gene_data.append(parent.gene_data[_])

            if i_check in switch_point:
                # 유전자를 받아오는 부모 변경
                try:
                    parent = parents[parents.index(parent) + 1]
                except IndexError:
                    parent = parents[0]

                """
                a = parents.index(parent) --> 현재 부모의 index 값
                parents[a+1] --> 부모 리스트에서, 현재 부모 인덱스값보다 +1 된 값 가져옴
                IndexError --> 만약 1이 넘어가면
                parent = parents[0] 다시 0으로 돌아옴
                """

        # return DNA(gene_data)
        dna = DNA(gene_data)
        return dna

    def evolution(self):
        print("Start Evolution Generation level %d" % Generation.cnt)

        dna_list = [self.best for _ in range(GOOD_DNA_CNT)]

        dna_list += [self.make_child() for _ in range(len(self.DNA_list) - len(dna_list))]

        return Generation(dna_list)

    @property
    def fitness(self):
        # 세대 객체의 평균 적합도
        return mean([dna.fitness for dna in self.DNA_list])

    @property
    def best(self):
        _best_dna = sorted(self.DNA_list, key=lambda x: x.fitness, reverse=True)[0]
        print(self.DNA_list[0])

        return _best_dna


class DNA:
    def __init__(self, gene_data=None):
        # 유전자 정보

        if gene_data is None:
            gen_list = []
            ran_num = randint(_gen_min_value, _gen_max_value)

            for i in range(_gen_len):
                while ran_num in gen_list:
                    ran_num = randint(_gen_min_value, _gen_max_value)
                gen_list.append(ran_num)

            self.gene_data = gen_list

        else:
            self.gene_data = gene_data

    def __repr__(self):
        return "< Gene %s | %d >" % ("_".join(str(x) for x in self.gene_data), self.fitness)

    @staticmethod
    def max_fitness():
        # 유전자 표기법 : 이진법
        _max_fitness = 100

        return _max_fitness

    @property
    def fitness(self) -> int:
        """
        적합도 계산 함수
        :return: 적합도 값
        """
        # max_jsgmd = (jsgmd_input + 1) * (jsgmd_input * 2) + (jsgmd_input * 2)
        # turtle = list(jsgmd_input * jsgmd_input)

        max_score = DNA.max_fitness()   # max_fitness : 100

        turtle_1 = self.gene_data[29] + self.gene_data[0] + self.gene_data[3] + self.gene_data[28] + self.gene_data[25] + self.gene_data[4]
        turtle_2 = self.gene_data[28] + self.gene_data[1] + self.gene_data[4] + self.gene_data[27] + self.gene_data[24] + self.gene_data[5]
        turtle_3 = self.gene_data[25] + self.gene_data[4] + self.gene_data[7] + self.gene_data[24] + self.gene_data[21] + self.gene_data[8]

        turtle_4 = self.gene_data[27] + self.gene_data[2] + self.gene_data[5] + self.gene_data[26] + self.gene_data[23] + self.gene_data[6]
        turtle_5 = self.gene_data[24] + self.gene_data[5] + self.gene_data[8] + self.gene_data[23] + self.gene_data[20] + self.gene_data[9]
        turtle_6 = self.gene_data[21] + self.gene_data[8] + self.gene_data[11] + self.gene_data[20] + self.gene_data[17] + self.gene_data[12]

        turtle_7 = self.gene_data[23] + self.gene_data[6] + self.gene_data[9] + self.gene_data[22] + self.gene_data[19] + self.gene_data[10]
        turtle_8 = self.gene_data[20] + self.gene_data[9] + self.gene_data[12] + self.gene_data[19] + self.gene_data[16] + self.gene_data[13]
        turtle_9 = self.gene_data[19] + self.gene_data[10] + self.gene_data[13] + self.gene_data[18] + self.gene_data[15] + self.gene_data[14]

        gap_1_2 = abs(turtle_1 - turtle_2)
        gap_1_3 = abs(turtle_1 - turtle_3)
        gap_1_4 = abs(turtle_1 - turtle_4)
        gap_1_5 = abs(turtle_1 - turtle_5)
        gap_1_6 = abs(turtle_1 - turtle_6)
        gap_1_7 = abs(turtle_1 - turtle_7)
        gap_1_8 = abs(turtle_1 - turtle_8)
        gap_1_9 = abs(turtle_1 - turtle_9)

        tmp_score = gap_1_2 + gap_1_3 + gap_1_4 + gap_1_5 + gap_1_6 + gap_1_7 + gap_1_8 + gap_1_9

        score = max_score - tmp_score

        return score


def visualization(generations):
    fitness_list = [generation.fitness for generation in generations]

    # 최대 적합도를 그래프에 나타냄
    max_fitness = DNA.max_fitness()
    plot([max_fitness for _ in range(len(generations))])

    xlim([0, len(generations)])

    # 축의 lim 값을 데이터 보다 높게 잡아줌으로써, 그래프의 가독성을 높임
    ylim([int(min(fitness_list)), (DNA.max_fitness() * 1.2)])

    xlabel('Generation')
    ylabel('Fitness Score')

    # 각 세대의 (평균) 적합도를 이용해 그래프에 나타냄
    plot([generation.fitness for generation in generations])

    show()


if __name__ == '__main__':
    ### input parameter ###
    gen_num = 100
    max_train = 200
    #######################

    Generations = list()

    # 첫 세대 (조상 세대)
    Generations.append(Generation([DNA() for _ in range(gen_num)]))

    i_check = 0
    for j in range(0, max_train):
        try:
            print("### episode : %d ###" % max_train)
            next_generation = Generations[i_check].evolution()
            Generations.append(next_generation)

            print("Fitness: %d" % next_generation.fitness)

            # 적합도가 최대일 경우, 반복문 종료
            if next_generation.fitness >= DNA.max_fitness() or j==(max_train-1):
                print("Best DNA: %s" % next_generation.best)
                break

            i_check += 1

        except KeyboardInterrupt:
            break
    
    print("Last Generation's Best DNA: %s" % Generations[-1].best)

    # print("---%d---" % Generations[-1].best[29])
    # print("--%d-%d--" % Generations[-1].best[0], Generations[-1].best[3])
    # print("--%d-%d--" % Generations[-1].best[28], Generations[-1].best[25])
    # print("-%d-%d-%d-" % Generations[-1].best[1], Generations[-1].best[4], Generations[-1].best[7])
    # print("-%d-%d-%d-" % Generations[-1].best[27], Generations[-1].best[24], Generations[-1].best[21])
    # print("%d-%d-%d-%d" % Generations[-1].best[2], Generations[-1].best[5], Generations[-1].best[8], Generations[-1].best[11])
    # print("%d-%d-%d-%d" % Generations[-1].best[26], Generations[-1].best[23], Generations[-1].best[20], Generations[-1].best[17])
    # print("-%d-%d-%d-" % Generations[-1].best[6], Generations[-1].best[9], Generations[-1].best[12])
    # print("-%d-%d-%d-" % Generations[-1].best[22], Generations[-1].best[19], Generations[-1].best[16])
    # print("--%d-%d--" % Generations[-1].best[10], Generations[-1].best[13])
    # print("--%d-%d--" % Generations[-1].best[18], Generations[-1].best[15])
    # print("---%d---" % Generations[-1].best[14])

    visualization(Generations)

