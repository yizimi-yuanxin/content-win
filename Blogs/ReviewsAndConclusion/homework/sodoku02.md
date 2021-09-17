### 〇、引言

大家久等的第二期来辣！

本期要讲述如何用程序来解答数独题目！实际上这个解法涉及算法，需要大家用一定的时间了解一下时间复杂度和搜索算法（与回溯）！

本期可能内容较复杂且枯燥，请耐心学习。

**本期关键词汇：暴力枚举算法、时间复杂度、可行性剪枝（广义）**

**本文为课后总结，除个人解释和思路外，内容均为上课老师讲解提供，请勿转载！！**

### 一、暴力枚举法枚举可能答案

众所周知电脑的运行速度很快（伪，理由见下文），所以对于一个数独题目，我们最暴力的方法就是直接把所有可能全部列出来，然后判断每种情况对不对。

我们先设置一个按钮ans，表示点击按钮就会自动填写正确答案：（具体讲解可参考上一期内容）

```py
ans = tk.Button(framebottom, text = "答案", relief = tk.RAISED,\
     font = ('HeiTi', '14', "bold"), width = 6, height = 1, bg = 'black', fg = 'white')
ans.pack(side = tk.LEFT, padx = 5)
# ===========================================
def ansclick(event):
    showlayout(solve(layout, fix), fix, gridvar)
ans.bind("<Button-1>", ansclick)
```
这里solve()是我们得出答案的函数。我们本文的重点就是写solve和其相关的函数。

首先，我们设置一个fix变量如下：
```py
fix = [[False for j in range(N)] for k in range(N)]
```
说明我们该行该列的数字是否是预先给出的数字，这样的话我们在枚举的时候保证它们不变。

我们假设一开始所有的空的答案都是$1$（必然不可能正确），然后我们找到第一个待填格子的坐标，逐个尝试和修改，从$1$到$9$。

```py
def solve(layout, fix):
    # 先把所有的待填数字改成1
    for row in range(N):
        for column in range(N):
            if (not fix[row][column]):
                layout[row][column] = 1
    # startgrid表示我们第一个待填空的位置，然后我们从这个空开始一个空一个空的尝试
    startgrid = [0, 0]
    for index in range(N * N):
        # 这里直接用81个格子统一编号来的，实际上可以上两重循环来找
        if not fix[index // N][index % N]:
            startgrid = [index // N, index % N]
            break
    
    solved = False
    # solved表示我们是否找到答案
    while (not solved) & (layout[startgrid[0]][startgrid[1]] <= 9):
        # 找到答案了自然就跳出循环，没找到答案且还没枚举完当下位置的答案
        # verify还是我们之前的判重函数（不过我们之后会进行优化，试想为什么(1)
        if verify():
            solved = True
            break
        else:
            layout = getnext(layout, fix)
            # 这里的getnext是将下一种情况列出的函数，返回layout的下一种情况(2)
    return layout
```

###### (1) verify的优化

我们之前的verify函数是这样写的：
```py
# def verify():
#     return verifyrow() & verifycolumn() & verifyblock()
```

对于Python语言来说，我们计算 '&' 的过程中（这里是指三个值的位运算与，可以理解为三个的返回值都是True式子才得出True，否则是False），要把 '&' 左右都先算出来才可以。那如果我已经得知verifyrow()是False，那么我根本不需要再执行下一个函数而可直接返回False，如果我计算的多了，可能会导致程序运行较慢。所以我们如此写：

```py
def verify(layout):
    if not verifyrow(layout):
        return False
    elif not verifycolumn(layout):
        return False
    elif not verifyblock(layout):
        return False
    else:
        return True
```

如果return了函数运行就结束了，所以我们可以这样避免无用功。我这里加了layout，在第一部分中可以去掉，而直接调用全局变量。

###### (2) getnext

getnext函数是用来寻找下一个可能的答案的函数。我们主要的思路就是拿过来上一个情况的layout，将这个表的最后一个可变化的数字找到，将其加一，如果加到9了，就让其归一，将上一个可变的变化的数字加一，如果仍然到9了，就以此类推。
代码如下：

```py
def getnext(layout, fix):
    nextchanged = False
    # 是否变化完毕
    for row in range(N - 1, -1, -1):
        # 这里最后一个-1是指自增的数量，-1就是倒序遍历
        for column in range(N - 1, -1, -1):
            if (not fix[row][column]):
                layout[row][column] += 1
                # 无论如何先自增一
                if layout[row][column] > 9:
                    layout[row][column] = 1
                    # 如需进位将其归一
                else :
                    nextchanged = True
                    # 如果没有出现进位或者进位完成，就直接退出
                    break
        if nextchanged:
            break
    return layout
```

### 二、时间复杂度

我们这就是用最简单最暴力的方法来枚举出所有可能的情况，然后判断其对错。当然你可以很清楚的认识到，这个方法狠傻大笨。原因很简单，我寻找了一大堆奇怪且不符合数独规则的数字，如果不考虑已经预先给出会减少的可能性，一共有多少种情况嘞？

答案是$9^{81}$次方。尽管我们电脑很快，但实际一秒仅可以运行计算约1亿次，假设我们一种情况只需要1次计算，我们需要计算完这些情况就大概需要$6.23 * 10^{61}$年的时间。这个时间显然太长。我们优化代码减少运行时间嘞？


在计算机科学中，算法的时间复杂度是一个函数，它定性描述该算法的运行时间。这是一个代表算法输入值的字符串的长度的函数。时间复杂度常用大$O$符号表述，不包括这个函数的低阶项和首项系数。使用这种方式时，时间复杂度可被称为是渐近的，亦即考察输入值大小趋近无穷时的情况。例如，如果一个算法对于任何大小为 $n$ （必须比 $n_0$ 大）的输入，它至多需要 $5n^3 + 3n$ 的时间运行完毕，那么它的渐近时间复杂度是 $O(n^3)$。[1]

在只有较大循环的程序中，我们一般考虑我们程序最大嵌套数的与$n$相关的多层嵌套循环作为时间复杂度的计算依据；如果我们有函数的递归，程序函数反复利用的话，就要多层考虑其计算了。

例如我们计算刚刚时间复杂度，我们一般设一行、一列的个数为$n = 9$，那我们考虑每一种情况的枚举，我们发现考虑一个位置的时间复杂度可能是$O(n)$，一共有$n^2$个格子，则其时间复杂度为$O(n^{n^{2}})$

之所以说这个时间复杂度巨大，是因为这个函数的递增速度奇快无比。$n = 1$时就是$1$，$n = 2$时就已经是$16$了，$n = 3$时就达到了$3^9$，以这个增长速度，我们电脑能够承受的运行顶多能运行$n = 3$的情况。

当然我们不可能用这样一种程序来计算数独，~~那样还不如找人做或者自己做~~，所以我们要另想方法。

### 三、挑选可行可能

我们其实可以从暴力做法进行优化。首先，我们可以先把该行该列预先给出的数字在每个位置的可能中删除（这里我们给每个位置一个列表，代表可能的数值，我们预先将其初始化为$[1,2,3,4,5,6,7,8,9]$)。

代码如下：
```py
def findpossiblevalues(layout, fix):
    values = copy.deepcopy(layout)
    # 这里就不能直接用等于了，否则原列表要被改变
    # 这里用深度copy
    for row in range(N):
        for column in range(N):
            # 枚举列表中的每一个位置
            if (not fix[row][column]):
                # 如果没有预先填出，则有1~9的可能
                values[row][column] = [i + 1 for i in range(N)]
            else:
                # 如果有预先填出，则可能性只有一个数字
                values[row][column] = [layout[row][column]]
    for row in range(N):
        for column in range(N):
            # 再次枚举列表位置，当为预填数字时
            if fix[row][column]:
                # print(values[row][column])
                number = values[row][column][0]
                # 将此数字挑出
                for i in range(N):
                    # 枚举该数字所在行的其他数字，在其可能中删除此数字
                    if (i != column) & (number in values[row][i]):
                        values[row][i].remove(number)
                for i in range(N):
                    # 枚举该数字所在列的其他数字，在其可能中删除此数字
                    if (i != row) & (number in values[i][column]):
                        values[i][column].remove(number)
                block = getblock_ij(row, column)
                # 寻找九宫的坐标，这里可以参考之前的找九宫的函数
                for i in range(block[0], block[1] + 1):
                    for j in range(block[2], block[3] + 1):
                        # 枚举该数字所在大九宫的其他数字，在其可能中删除此数字
                        if ((i != row) | (j != column)) & (number in values[i][j]):
                            values[i][j].remove(number)
    return values
```

这时理论的时间复杂度还是很大，但是实际可以运行的$n$会变大很多，试想为什么。

这时，我们solve函数和getnext函数也会有所变化。我们本来solve函数是先把没有预先填出的位置初始化成1，我们现在可以初始化为第一个可能的值；相同的，getnext函数本来是将其自增1，我们可以改为将其寻找下一个可能值，或者归为第一个可能值。代码如下：（仅注释修改部分的代码）

```py
def getnext(layout, fix, values):
    nextchanged = False
    for row in range(N - 1, -1, -1):
        for column in range(N - 1, -1, -1):
            if (not fix[row][column]):
                currentvalues = values[row][column]
                # 有正确的可能的列表传递
                i = currentvalues.index(layout[row][column])
                # 寻找当前layout最后位置的数字在可能的列表中的位置
                if i < (len(currentvalues) - 1): # 如果在正确的位置之中
                    layout[row][column] = currentvalues[i + 1]
                    nextchanged = True
                    break
                else: # 如果位置信息超过正常值，就是枚举完了，就重新进行类似归一的操作
                    layout[row][column] = currentvalues[0]
        if nextchanged:
            break
    return nextchanged, layout

def solve(layout, fix):
    values = findpossiblevalues(layout, fix)
    # 先找出一个可能值的列表，注意这里是三维，
    # 第一维是行，第二维是列，第三维是该行该列可能值，
    # 所以尽管可能值只有一个，也不能直接调用两维而忽略其为三维列表
    for row in range(N):
        for column in range(N):
            if not fix[row][column]:
                layout[row][column] = values[row][column][0]
                # 初始化为第一个可能值
    
    solved = False
    nextchanged = True
    while (not solved) & (nextchanged):
        if verify(layout):
            # 通过则给出答案
            solved = True
            break
        else:
            nextchanged, layout = getnext(layout, fix, values)
            # 不通过则寻找下一个值，当然如果没有下一个值也要跳出
    return layout
    # 返回答案
```

###### 总体代码一览：
```py
'''
writer : yizimi - yuanxin
Instructor : Mr. Mao, Place of Tang Dynasty and CITers
'''

import tkinter as tk
import numpy as np
from tkinter.messagebox import showinfo
import copy

N = 9
layout = [[0 for j in range(N)] for k in range(N)]
# ======================================================
fix = [[False for j in range(N)] for k in range(N)]
# ======================================================

root = tk.Tk()
root.title("数独游戏")

frametop = tk.Frame(root)
gridvar = [[tk.StringVar() for column in range(N)] for row in range(N)]
frame = [tk.Frame(frametop) for row in range(N)]
grid = [[tk.Button(frame[row], width = 3, textvariable = gridvar[row][column],\
     relief = tk.GROOVE, command = lambda row = row,\
          column = column:gridclick(row, column),\
          font = ('Helvetica', '12')) for column in range(N)] for row in range(N)]
for row in range(N):
    for column in range(N):
        grid[row][column].pack(side = tk.LEFT)
    frame[row].pack(side = tk.TOP)
frametop.pack(side=tk.TOP, pady = 10)

framemiddle = tk.Frame(root)
selections = [tk.Button(framemiddle, width = 3, text = '%d' % number, relief = tk.RAISED,\
     command = lambda number=number:numberclick(selections[number - 1]),\
          font = ('Helvetica', '12')) for number in range(1, 10)]
for each in selections:
    each.pack(side = tk.LEFT)
framemiddle.pack(side=tk.TOP, pady = 15)


framebottom = tk.Frame(root)
erase = tk.Button(framebottom, text = '删除', relief = tk.RAISED,\
     font = ('HeiTi', '14', 'bold'), width = 6, height = 1, bg = 'darkgreen', fg = 'white')
erase.pack(side = tk.LEFT, padx = 5)
check = tk.Button(framebottom, text = '核查', relief = tk.RAISED,\
     font = ('HeiTi', '14', 'bold'), width = 6, height = 1, bg = 'darkblue', fg = 'white')
check.pack(side = tk.LEFT, padx = 5)
ans = tk.Button(framebottom, text = "答案", relief = tk.RAISED,\
     font = ('HeiTi', '14', "bold"), width = 6, height = 1, bg = 'black', fg = 'white')
ans.pack(side = tk.LEFT, padx = 5)
ok = tk.Button(framebottom, text = '退出', relief = tk.RAISED, command = exit,\
     font = ('HeiTi', '14', 'bold'), width = 6, height = 1, bg = 'darkred', fg = 'white')
ok.pack(side = tk.LEFT, padx = 5)
framebottom.pack(side = tk.TOP, pady = 5)

def gridclick(row, column):
    number = ''
    for i in range(N):
        if selections[i]["relief"] == tk.SOLID:
            number = '%d' % (i + 1)
            break
    gridvar[row][column].set(number)
    if number == '':
        layout[row][column] = 0
    else :
        layout[row][column] = int(number)


def numberclick(selectionbutton):
    for i in range(N):
        selections[i]['relief'] = tk.RAISED
    erase["relief"] = tk.RAISED
    selectionbutton["relief"] = tk.SOLID

def eraseclick(event):
    for i in range(N):
        selections[i]['relief'] = tk.RAISED


def checkclick(event):
    correct = verify(layout) 
    if correct:
        showinfo("核查结果", "答案正确")
        print("correct")
    else:
        showinfo("核查结果", "答案不正确")
        print("wrong")
    # check['relief'] = tk.RAISED
check.bind("<Button-1>", checkclick)

# ========================================================
def ansclick(event):
    showlayout(solve(layout, fix), fix, gridvar)
ans.bind("<Button-1>", ansclick)
# ========================================================

def readlayout(filename):
    layoutfile = open(filename, 'r')
    lines = layoutfile.readlines()
    for row in range(N):
        line = lines[row].strip('\n')
        if line != '':
            for i in range(len(line)):
                if line[i] != '' and line[i] != ' ':
                    layout[row][i] = int(line[i])
    return layout

def showlayout(layout, fix, gridvar):
    for row in range(N):
        for column in range(N):
            gridvar[row][column].set(str(layout[row][column]) if (layout[row][column] != 0) else '')
            if fix[row][column]:
                grid[row][column]["state"] = tk.DISABLED


def verifyrow(layout):
    correct = True
    for row in range(N):
        line = layout[row].copy()
        # print(line)
        line.sort()
        if line != [1, 2, 3, 4, 5, 6, 7, 8, 9]:
            correct = False
    return correct

def verifycolumn(layout):
    correct = True
    for column in range(N):
        line = list((np.array(layout))[:, column])
        line.sort()
        if line != [1, 2, 3, 4, 5, 6, 7, 8, 9]:
            correct = False
    return correct

def verifyblock(layout):
    correct = True
    for blockindex in range(N):
        block = getblock(blockindex)
        line = list((np.array(layout))[block[0] : block[1] + 1, block[2] : block[3] + 1].reshape(N))
        line.sort()
        if line != [1, 2, 3, 4, 5, 6, 7, 8, 9]:
            correct = False
    return correct


def getblock(index):
    rowstart = index // 3 * 3
    rowend = rowstart + 2
    columnstart = index % 3 * 3
    columnend = columnstart + 2
    return rowstart, rowend, columnstart, columnend

# =====================必要性添加========================

# def verify():
#     return verifyrow() & verifycolumn() & verifyblock()
def verify(layout):
    if not verifyrow(layout):
        return False
    elif not verifycolumn(layout):
        return False
    elif not verifyblock(layout):
        return False
    else:
        return True

def findfix(layout):
    fix = [[False for j in range(N)] for k in range(N)]
    for row in range(N):
        for column in range(N):
            if layout[row][column] > 0:
                fix[row][column] = True
    return fix

# =======================================================

# =======================暴力枚举法======================= 

# def getnext(layout, fix):
#     nextchanged = False
#     for row in range(N - 1, -1, -1):
#         for column in range(N - 1, -1, -1):
#             if (not fix[row][column]):
#                 layout[row][column] += 1
#                 if layout[row][column] > 9:
#                     layout[row][column] = 1
#                 else :
#                     nextchanged = True
#                     break
#         if nextchanged:
#             break
#     return layout

# def solve(layout, fix):
#     for row in range(N):
#         for column in range(N):
#             if (not fix[row][column]):
#                 layout[row][column] = 1
#     startgrid = [0, 0]
#     for index in range(N * N):
#         if not fix[index // N][index % N]:
#             startgrid = [index // N, index % N]
#             break
    
#     solved = False
#     while (not solved) & (layout[startgrid[0]][startgrid[1]] <= 9):
#         if verify():
#             solved = True
#             break
#         else:
#             layout = getnext(layout, fix)
#     return layout

# ==========================================================

# ========================枚举剪枝法=========================

def findpossiblevalues(layout, fix):
    values = copy.deepcopy(layout)
    for row in range(N):
        for column in range(N):
            if (not fix[row][column]):
                values[row][column] = [i + 1 for i in range(N)]
            else:
                values[row][column] = [layout[row][column]]
    for row in range(N):
        for column in range(N):
            if fix[row][column]:
                # print(values[row][column])
                number = values[row][column][0]
                for i in range(N):
                    if (i != column) & (number in values[row][i]):
                        values[row][i].remove(number)
                for i in range(N):
                    if (i != row) & (number in values[i][column]):
                        values[i][column].remove(number)
                block = getblock_ij(row, column)
                for i in range(block[0], block[1] + 1):
                    for j in range(block[2], block[3] + 1):
                        if ((i != row) | (j != column)) & (number in values[i][j]):
                            values[i][j].remove(number)
    return values

def getblock_ij(row, column):
    rowstart = row // 3 * 3
    rowend = rowstart + 2
    columnstart = column // 3 * 3
    columnend = columnstart + 2
    return rowstart, rowend, columnstart, columnend

def getnext(layout, fix, values):
    nextchanged = False
    for row in range(N - 1, -1, -1):
        for column in range(N - 1, -1, -1):
            if (not fix[row][column]):
                currentvalues = values[row][column]
                i = currentvalues.index(layout[row][column])
                if i < (len(currentvalues) - 1):
                    layout[row][column] = currentvalues[i + 1]
                    nextchanged = True
                    break
                else:
                    layout[row][column] = currentvalues[0]
        if nextchanged:
            break
    return nextchanged, layout

def solve(layout, fix):
    values = findpossiblevalues(layout, fix)

    for row in range(N):
        for column in range(N):
            if not fix[row][column]:
                layout[row][column] = values[row][column][0]
    
    solved = False
    nextchanged = True
    while (not solved) & (nextchanged):
        if verify(layout):
            solved = True
            break
        else:
            nextchanged, layout = getnext(layout, fix, values)

    return layout

# ==========================================================


erase.bind("<Button-1>", eraseclick)

layout = readlayout('solution_revised.txt')
fix = findfix(layout)

showlayout(layout, fix, gridvar)

tk.mainloop()
```

### 附：后话

这个版本仍然是很慢的方法，我们可以从判断可能的情况中数字重复或不符合规则的可能来把更多错误情况加速程序。这个涉及**搜索算法和回溯算法**，推荐在下期更新之前预习一下！QwQ


参考文献及网页：
    [0]老师上课的精彩讲解及老师的资料
    [1]https://baike.sogou.com/v105680.htm?fromTitle=%E6%97%B6%E9%97%B4%E5%A4%8D%E6%9D%82%E5%BA%A6
