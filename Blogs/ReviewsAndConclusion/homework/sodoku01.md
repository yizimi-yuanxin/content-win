### 〇、引言
QwQ

我们即将用 Python 写一个GUI图形界面数独！（第一部分）

设计效果：
![hjtiy-vc7oi.gif](https://i.loli.net/2020/12/10/hxg8JtRAzODqTBZ.gif)

**关键词汇：tkinter库、python方法判重、GUI界面简单设计**

**本文为课后总结，除个人解释和思路外，内容均为上课老师讲解提供，请勿转载！！**

### 一、tkinter库及数独界面设计

##### 1.GUI界面创建

tkinter，一个神奇的东西。Python自带的控件，只需要调用：

```py
import tkinter as tk
```

我们首先要创建一个图形界面的根界面。我们可以：
```py
root = tk.Tk()
root.title("数独游戏")
```

此时root是一个GUI界面的类，root.title()就是给图形界面加标题。

效果：无！

当你创建这个界面的时候，你会发现这个揭=界面很快就会被关掉（电脑好一点的话根本看不到它出现），这个时候你需要让这个界面保持工作，也就是我们让他进入服务器式的循环中不关闭。我们可以将这行代码加在后面：

```py
tk.mainloop()
```

此时效果：
![](https://cdn.luogu.com.cn/upload/image_hosting/9doxbf3b.png)

~~非常朴素的开始~~

##### 2.初识控件

你看到的一些奇怪的文字框，按钮啥的，都是用控件组成。我们所看到上文的演示的控件主要是按钮控件（Button）。我们也相应的介绍一下其他常用的控件：

###### (1) label

只是一个文字框

```py
label = tk.Label(root, fg='red', bg='blue',\
     width=10, height=2, text='标签示例', font=('Tempus Sans ITC', 12))
label.pack()
```

    fg: 字体颜色
    bg: 背景颜色
    width: 宽度
    height: 高度
    text: 文字内容
    font: 字体（字体，字号）

label属于tk.Label类，需要用label.pack()打包然后在界面中显示。显示情况如下：
![](https://cdn.luogu.com.cn/upload/image_hosting/qykcxeaz.png?x-oss-process=image/resize,m_lfit,h_170,w_225)

###### (2) entry

可输入文字框

```py
pwd = tk.StringVar()
entry = tk.Entry(root, textvariable=pwd, relief=tk.RAISED)
pwd.set("输入框示例")
entry.pack()
```

    textvariable: 框中初始所填的文字
    relief: 形态/状态

在使用界面编程的时候，有些时候是需要跟踪变量的值的变化，以保证值的变更随时可以显示在界面上。由于python无法做到这一点，所以使用了tcl的相应的对象，也就是StringVar、BooleanVar、DoubleVar、IntVar所需要起到的作用[1]。我们这里用的是StringVar()。

relief表示形态或者状态，其实就是图形框的样式。常见的有"FLAT", "RAISED", "SUNKEN", "SOLID", "RIDGE", "GROOVE"，但注意修改时这些变量都是tkinter内部的。其样式效果如下（下面是用的按钮展示的）：

![](https://cdn.luogu.com.cn/upload/image_hosting/ctuv6jiv.png)

Entry控件的演示：
![Entrykongjian.gif](https://i.loli.net/2020/12/10/H6ODQ89dMJliYZC.gif)

###### (3) Botton

按钮，可用于执行命令

```py
def buttonclicked():
    return True
btn = tk.Button(root, text="按钮示例", relief=tk.SOLID, bd=2, command=buttonclicked).pack()
```

    text: 显示的文字
    relief: 形态/样式
    bd: 按钮的边缘宽度(borderwidth)
    command: 回调函数，当你按下这个按钮后会执行的函数
    font: 文字的字体和字号
    width: 宽度
    height: 高度
    bg: 背景颜色
    fg: 字体颜色

效果:
![](https://cdn.luogu.com.cn/upload/image_hosting/kxdtbim8.png)

代码：
```py
import tkinter as tk

root = tk.Tk()
root.title("数独游戏")

#Label
label = tk.Label(root, fg='red', bg='blue', width=10, height=2, text='标签示例', font=('Tempus Sans ITC', 12))
label.pack()

#Entry
pwd = tk.StringVar()
entry = tk.Entry(root, textvariable=pwd, relief=tk.RAISED)
pwd.set("输入框示例")
entry.pack()

#Button
def buttonclicked():
    return True
btn = tk.Button(root, text="按钮示例", relief=tk.SOLID, bd=2, command=buttonclicked).pack()

tk.mainloop()
```

###### (4) 更多

更多控件：

    Button	        按钮控件；在程序中显示按钮。
    Label	        标签控件；可以显示文本和位图
    Entry	        输入控件；用于显示简单的文本内容
    Text	        文本控件；用于显示多行文本
    Radiobutton	单选按钮控件；显示一个单选的按钮状态
    Checkbutton	多选框控件；用于在程序中提供多项选择框
    Listbox	        列表框控件；在Listbox窗口小部件是用来显示一个字符串列表给用户
    Frame	        框架控件；在屏幕上显示一个矩形区域，多用来作为容器
    Canvas	        画布控件；显示图形元素如线条或文本
    Menubutton	菜单按钮控件，由于显示菜单项
    Menu	        菜单控件；显示菜单栏,下拉菜单和弹出菜单
    Message	        消息控件；用来显示多行文本，与label比较类似
    Scale	        范围控件；显示一个数值刻度，为输出限定范围的数字区间
    Scrollbar	滚动条控件，当内容超过可视化区域时使用，如列表框
    Toplevel	容器控件；用来提供一个单独的对话框，和Frame比较类似
    Spinbox	        输入控件；与Entry类似，但是可以指定输入范围值
    PanedWindow	PanedWindow是一个窗口布局管理的插件，可以包含一个或者多个子控件
    LabelFrame	labelframe 是一个简单的容器控件。常用于复杂的窗口布局
    tkMessageBox	用于显示你应用程序的消息框

更多标准属性：

    属性	    描述
    Dimension   控件大小
    Color	    控件颜色
    Font	    控件字体
    Anchor	    锚点
    Relief	    控件样式
    Bitmap	    位图
    Cursor	    光标

更多集合管理：

    几何方法    描述
    pack()	    包装
    grid()	    网格
    place()	    位置

##### 3.数独界面设计

我们发现，我们目标效果的界面可以大体分为三个界面：

![](https://cdn.luogu.com.cn/upload/image_hosting/n3rf03ck.png)

对于大的分区我们可以用Frame控件来调整和分区。

Python Tkinter 框架（Frame）控件在屏幕上显示一个矩形区域，多用来作为容器[2]。我们可以用Frame框架来分区。这时，我们根据刚刚的分区效果来看，我们可以这样写：

```py
import tkinter as tk

root = tk.Tk()
root.title("数独游戏")

frametop = tk.Frame(root)
# 上部框架建设
frametop.pack(side=tk.TOP, pady = 10)
# side指该框架要放在哪里，我们可以选择TOP, BOTTOM, LEFT, RIGHT
# pady是指与垂直边距，相似的，padx是指水平边距

framemiddle = tk.Frame(root)
framemiddle.pack(side=tk.TOP, pady=15)

framebottom = tk.Frame(root)
framebottom.pack(side=tk.TOP, pady=5)

tk.mainloop()
```

效果：空界面

因为此时我们并没有在各个框架上添加内容。我们可以通过加控件来丰富我们的框架。对于我们数独游戏界面来说，我们可以这样写：(部分代码解释见备注，控件表示请见上文表格)

```py
import tkinter as tk

root = tk.Tk()
root.title("数独游戏")

N = 9

frametop = tk.Frame(root)
# 上部的框架建设
gridvar = [[tk.StringVar() for column in range(N)] for row in range(N)]
# 对于大九宫格的所有的显示都是大多动态的，我们分别定义一个StringVar()来存储
frame = [tk.Frame(frametop) for row in range(N)]
# 这里涉及框架的嵌套，请见下文详解(1)
grid = [[tk.Button(frame[row], width = 3, textvariable = gridvar[row][column],\
     relief = tk.GROOVE, command = lambda row = row,\
          column = column:gridclick(row, column),\
          font = ('Helvetica', '12')) for column in range(N)] for row in range(N)]
# grid这里为控件列表，代表每个小格对应的Button控件，lambda的解释请看下文详解(2)
for row in range(N):
    for column in range(N):
        grid[row][column].pack(side = tk.LEFT)
        # 分别对每个控件进行打包显示
    frame[row].pack(side = tk.TOP)
    # 对嵌套的框架进行打包显示
frametop.pack(side=tk.TOP, pady = 10)
# 大框架进行打包显示

framemiddle = tk.Frame(root)
# 中部的框架建设
selections = [tk.Button(framemiddle, width = 3, text = '%d' % number, relief = tk.RAISED,\
     command = lambda number=number:numberclick(selections[number - 1]),\
          font = ('Helvetica', '12')) for number in range(1, 10)]
# 设立九个“选项”按钮
for each in selections:
    each.pack(side = tk.LEFT)
framemiddle.pack(side=tk.TOP, pady = 15)
# 层层打包显示

framebottom = tk.Frame(root)
# 底部的框架建设
erase = tk.Button(framebottom, text = '删除', relief = tk.RAISED,\
     font = ('HeiTi', '14', 'bold'), width = 7, height = 1, bg = 'darkgreen', fg = 'white')
# 删除键
erase.pack(side = tk.LEFT, padx = 15)
check = tk.Button(framebottom, text = '核查', relief = tk.RAISED,\
     font = ('HeiTi', '14', 'bold'), width = 7, height = 1, bg = 'darkblue', fg = 'white')
# 核查键
check.pack(side = tk.LEFT, padx = 15)
ok = tk.Button(framebottom, text = '退出', relief = tk.RAISED, command = exit,\
     font = ('HeiTi', '14', 'bold'), width = 7, height = 1, bg = 'darkred', fg = 'white')
# 退出键
ok.pack(side = tk.LEFT, padx = 15)
framebottom.pack(side = tk.TOP, pady = 5)
# 依次打包

tk.mainloop()
```

详解：

(1)
```py
frametop = tk.Frame(root)
frame = [tk.Frame(frametop) for row in range(N)]
```
这两行代码第一行是在root中添加一个框架，而第二行是在root下的一个框架中添加一个子框架，我们可以将其视作嵌套框架。这样方便我们对button们进行排版。框架建构也可以进行批量操作。

(2)
```py
grid = [[tk.Button(frame[row], width = 3, textvariable = gridvar[row][column],\
     relief = tk.GROOVE, command = lambda row = row,\
          column = column:gridclick(row, column),\
          font = ('Helvetica', '12')) for column in range(N)] for row in range(N)]
```
lambda：lambda是定义了匿名函数。一般我们用lambda来定义单行函数比如

```py
>>> lambda x : x + 1 (1)
2
```
这里第一个x代表函数的变量，冒号后面表示函数表达式或者单行函数所要调用的东西，也就是说这里的用法相当于定义了：
```py
def g(x):
    return x + 1 
```
我们源代码上如此用法实际上是为了方便调用函数，可以用方括号内所定义的变量来作为形参调用函数。

不过在上述代码中，我们numberclick，gridclick函数还未编写，按钮的功能函数还未完善。~~革命尚未完成，同志还需努力~~

### 二、内部构造架构

##### 1.numberclick -- 选项架构

我们的期望效果是点击选项架构后，再点击大九宫格中的某个格子时，将选项中的数字填入其大九宫格格子中。那我们可以让程序记住我们选项中选的数字，然后再填入。我们numberclick就是用来记录选中的数字的

```py
def numberclick(selectionbutton):
    # selectionbutton是我们点中的Button类
    for i in range(N):
        selections[i]["relief"] = tk.RAISED
    # 先将所有数字按钮和控制按钮恢复状态
    erase["relief"] = tk.RAISED
    # 恢复删除按钮显示状态
    selectionbutton["relief"] = tk.SOLID
    # 将点击的数字按钮设置成SOLID
```
这时，我们就将我们选中的格子用SOLID的方式记录下来，方便我们填数

##### 2.gridclick -- 大九宫格架构

选择所需填的数以后，我们在点击大九宫格的格子时，就可以将标记过的选项数字填入其中。

```py
def gridclick(row, column):
    number = ''
    # 一般首次点击选项之前都没有可选的number，那么初始化为''
    for i in range(N):
        if selections[i]["relief"] == tk.SOLID:
            number = '%d' % (i + 1)
            break
    # 这一for循环寻找选项中的标记项，然后将其下标+1（下标为0-8，我们数字实为1-9）作为待填项，注意我们之前用的时StringVar()，所有我们此时要修改也是str类型的，所以转换成str
    gridvar[row][column].set(number)
    # 将大九宫格对应行列的StringVar()类改为number
    if number == '':
        layout[row][column] = 0
    # 当然，如果没有数的话强制转换是肯定不行的
    else :
        layout[row][column] = int(number)
    # 这里layout表示现在的情况，这是我们定义的全局变量，以int记录，方便计算和判重
```
我们把上面两个函数放入代码中，运行效果如下：

![gongneng1.gif](https://i.loli.net/2020/12/10/zRSlDNdIb7pnjXM.gif)

##### 3.eraseclick -- 删除键功能架构

```py
def eraseclick(event):
    for i in range(N):
        selections[i]["relief"] = tk.RAISED
    # 点击删除后会将所有标记去掉，这样在填数的时候number变量为空
erase.bind("<Button-1>", eraseclick)
```

erase 是我们上文的所讲的Button控件，bind是将其他函数和控件绑定在一起的函数。其参数中 "\<Button\-1>" 表示左键点击， "\<Button\-2>" 指右键点击。也就是说，当我左键点击删除键时，程序运行eraseclick函数。

##### 4.checkclick -- 核查功能架构

```py
def checkclick(event):
    correct = verify()
    if correct:
        showinfo('核查结果', '答案正确') 
        print("答案正确")
    else:
        showinfo('核查结果', '答案不正确')  
        print("答案不正确")
check.bind("<Button-1>", checkclick)
```
注意：这里的最后一句是函数外的。这句话相当于初始化，一定不要将此句错缩进进函数中

verify()是一个返回Bool值的自定义函数，表示我们的填入是否是正确的。

这里涉及showinfo()函数，这是跳出一个提示窗口，其中形参第一个是**窗口名称**，第二个是**提示内容**，效果如下：
![showinfo_.gif](https://i.loli.net/2020/12/10/UdrpGcz2m4hTQM6.gif)

##### 4.readlayout -- 读入架构

我们做数独肯定不是一张空空的表格来让我们填的，而是有初始定下的几个数字。我们可以将题目提前存在文件中（这里我们将文件命名为"sodoku.txt"，储存方式如下：
```
8
  36     
 7  9 2  
 5   7   
    457
   1   3 
  1    68
  85   1
 9    4  
```

首先我们打开文件，这里我们直接将文件名储存在filename变量中作为参数。并将其以行读入成列表：
```py
    layoutfile = open(filename, 'r')
    lines = layoutfile.readlines()
```
随后我们逐行操作，先把每行的'\n'去掉，然后对行内每个字符进行统计和存储。注：可能存在一行中没有数字的可能，所以要看本行是否有数字，最后返回其数组。完整代码：

```py
def readlayout(filename):
    layoutfile = open(filename, 'r')
    lines = layoutfile.readlines()
    for row in range(N):
        line = lines[row].strip('\n')
        if line != '': # 除去空行情况
            for i in range(len(line)):
                if line[i] != '' and line[i] != ' ':
                    layout[row][i] = int(line[i])
    return layout
```

我们显示数字时，要把预先给出的数字做处理，使得其不能被改动。我们对每行每列的layout进行判定，若其中有数字，则将其性质该为不可变性（tk.DISABLED），其代码如下：
```py
def showlayout(layout, gridvar):
    for row in range(N):
        for column in range(N):
            gridvar[row][column].set(str(layout[row][column]) if (layout[row][column] != 0) else '') # 将已经有的预填入表中
            if layout[row][column] != 0:
                grid[row][column]["state"] = tk.DISABLED
```

### 三、核查答案正确性：判重

我们现在可以填数字了，现在我们需要对答案进行正确性核查。我们要保证在同行，同列和同九宫中不重复数字为1~9。首先我们要分别取到各行，各列，各九宫的数字。我们分别用3个函数来判定行，列与和九宫中数的正确性，其返回值为一个Bool变量。

我们每行用切片的方式切出，然后将9个数sort一下，排序后应该一定是1，2，3，4，5，6，7，8，9。

```py
def verifyrow(): # 按行取
    correct = True
    for row in range(N):
        line = layout[row].copy()
        line.sort()
        if line != [1, 2, 3, 4, 5, 6, 7, 8, 9]:
            correct = False
    return correct

def verifycolumn(): # 按列取
    correct = True
    for column in range(N):
        line = list((np.array(layout))[:, column])
        line.sort()
        if line != [1, 2, 3, 4, 5, 6, 7, 8, 9]:
            correct = False
    return correct

def verifyblock(): # 按九宫取
    correct = True
    for blockindex in range(N):
        block = getblock(blockindex)
        line = list((np.array(layout))[block[0] : block[1] + 1, block[2] : block[3] + 1].reshape(N))
        # 这里是切片，切出对应行和列。然后将其重塑成一个一维数组，再转成list方便sort
        line.sort()
        if line != [1, 2, 3, 4, 5, 6, 7, 8, 9]:
            correct = False
    return correct


def getblock(index): # 这个函数是用来获得区块所在的行的开始和结束，列的开始和结束
    rowstart = index // 3 * 3
    rowend = rowstart + 2
    columnstart = index % 3 * 3
    columnend = columnstart + 2
    return rowstart, rowend, columnstart, columnend
```

最终我们将行，列，九宫所得的正确性综合一下，确定最终核查结果，即：

```py
def verify():
    return verifyrow() & verifycolumn() & verifyblock()
```

这样，我们就大体将主要的效果搞出来了。所有代码综合起来如下：
```py
'''
writer : yizimi - yuanxin
Instructor : Mr. Mao, Place of Tang Dynasty and CITers
'''

import tkinter as tk
import numpy as np
from tkinter.messagebox import showinfo

N = 9
layout = [[0 for j in range(N)] for k in range(N)]

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
     font = ('HeiTi', '14', 'bold'), width = 7, height = 1, bg = 'darkgreen', fg = 'white')
erase.pack(side = tk.LEFT, padx = 15)
check = tk.Button(framebottom, text = '核查', relief = tk.RAISED,\
     font = ('HeiTi', '14', 'bold'), width = 7, height = 1, bg = 'darkblue', fg = 'white')
check.pack(side = tk.LEFT, padx = 15)
ok = tk.Button(framebottom, text = '退出', relief = tk.RAISED, command = exit,\
     font = ('HeiTi', '14', 'bold'), width = 7, height = 1, bg = 'darkred', fg = 'white')
ok.pack(side = tk.LEFT, padx = 15)
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
    correct = verify() 
    if correct:
        showinfo("核查结果", "答案正确")
        print("correct")
    else:
        showinfo("核查结果", "答案不正确")
        print("wrong")
check.bind("<Button-1>", checkclick)

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

def showlayout(layout, gridvar):
    for row in range(N):
        for column in range(N):
            gridvar[row][column].set(str(layout[row][column]) if (layout[row][column] != 0) else '')
            if layout[row][column] != 0:
                grid[row][column]["state"] = tk.DISABLED


def verifyrow():
    correct = True
    for row in range(N):
        line = layout[row].copy()
        line.sort()
        if line != [1, 2, 3, 4, 5, 6, 7, 8, 9]:
            correct = False
    return correct

def verifycolumn():
    correct = True
    for column in range(N):
        line = list((np.array(layout))[:, column])
        line.sort()
        if line != [1, 2, 3, 4, 5, 6, 7, 8, 9]:
            correct = False
    return correct

def verifyblock():
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

def verify():
    return verifyrow() & verifycolumn() & verifyblock()

erase.bind("<Button-1>", eraseclick)

layout = readlayout('sudoku.txt')
showlayout(layout, gridvar)

tk.mainloop()
```

##### # 但是还没有结束哦，我们下一期讲解如何自动填写正确答案QwQ

参考文献及网站：
    [0].老师精彩的课上讲解和资料（不方便透露其相关信息）
    [1].https://blog.csdn.net/Eider1998/article/details/104725180/
    [2].https://www.runoob.com/python/python-tk-frame.html