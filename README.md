#Integer Programming solver

Para executar: python integer_programming.py <<arquivo_de_entrada>>

O programa lê a entrada, trasnforma a PL em FPI e multiplica por -1 qualquer linha que tenha o b negativo.
Após isso, a PL auxiliar é criada e resolvida para verificar se a PL é viável. Caso seja e a gente tenha multiplicado
alguma linha da PL original por -1, o tableau da PL auxiliar é mantida para termos uma base viável pra inicar o simplex
na PL original.

O simplex primal então é executado na PL original e é verificado se há algum valor de x fracionário. Caso haja, chamamos
o plano de cortes ou o branch and bound, dependendo do parâmetro de entrada.

O plano de cortes foi implementado da seguinte maneira:
- em um loop "infinito" é verificado se há algum valor de x do problema original que é fracionário, é selecionada a primeira
- caso não haja, o ótimo e o valor de x da última execução do simplex dual é retornado
- caso haja, a linha correspondente a esse valor de x é selecionada para entrar como restrição
- é feito um floor em todas as entradas da linha e ela é adicionada no fim do tableau
- além disso é adicionada uma nova variável de folga para essa linha
- simplex dual é chamado para esse novo tableau
- repete o processo

O branch and bound foi implementado da seguinte maneira:
- é selecionado o primeiro de x do problema original que é fracionário (não existe a possibilidade de não ter um x que não é inteiro pela estrutura do programa)
- são criadas duas novas pl's (uma da "esquerda" e outra "direita") que são cópias da pl atual
- são adicionadas as restrições correspondentes (uma com floor de b e outra com ceil e com folga = -1)
- há uma variável global para guardar o melhor valor obj inteiro e o seu x correspondente
- o simplex dual é chamado para a pl esquerda
- se ela for viável e seu valor obj for maior que o da variável global, há duas opções:
	- se o seu x for inteiro para os x do problema original, então ela é a melhor inteira até o momento e as var globais são atualizadas
	- senão, é chamado o branch and bound pra essa nova pl recursivamente
- logo após, o processo é feito para a pl da direita
- portanto, fazemos uma busca em profundidade no espaço de soluções
