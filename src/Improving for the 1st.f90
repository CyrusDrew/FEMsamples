! this code is a open-source program for the simulation of a concentrate 
! force at end of a cantiever beam 
! Copyright (C) Yunrui Zhu 

! This code is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
! 
! You should have received a copy of the GNU General Public License
! along with this code.  If not, see <http://www.gnu.org/licenses/>.
! 											 		
! 											   |F 					
! \|___________________________________________↓
! \|	 								  	   |			
! \|		   triangle meshs				   | 	
! \|___________________________________________|
! \|
!
!	
!	Author:		Yunrui Zhu
!	Email:		yunruizhu@126.com
!	Version:	1.0
!	Created:	2014.05.02	Yunrui Zhu
!	Updated:	2015.07.11	Yunrui Zhu

program Improving4First

implicit none
    real,parameter :: e0=2.1e6,p0=0.3,Number_Nodes=85
    !e0，材料弹性模量；p0，材料泊松比
    integer :: NN,NE,NF,i,j,a,t,x
    !NN:节点总数，NE：单元总数，NF：方程总数（2*NN）
    real(8) :: B(3,6),D(3,3),K_e(6,6),inv_K_e(6,6),area
    !B：直角在左下角的三角形单元的B矩阵，invB是相反单元的B矩阵,
    !inv_K_e是反三角单元的刚度，理论上等于K_e
    !D:D矩阵
    integer :: E_N(128),E_i(128),E_j(128),E_m(128),Element_Num(128,3),El_N(3,128)
    !E_N:单元编号,E_i：单元的i节点（j,m同）,E_j,E_m;
    real(8) :: K(2*Number_Nodes,2*Number_Nodes)=0 !总体刚度矩阵
    integer :: G(6,170),G_Num,G_Element(2,2)
    !G:G转换矩阵；G_Num:用于计数；G_Element:G矩阵的单位阵
    real(8) :: wy(170),p(170),temp ,N_w(2,85)! wy位移，p载荷,N_w,节点位移
    real(8) :: Total_N_w(6,128),strain(3,128),stress(3,128),max_stress,min_stress
    !Total_N_w总的节点的位移集成矩阵
    real(8) :: e_a(2,2),e_b(2,2),e_c(2,2),e_d(2,2),e_e(2,2),e_f(2,2),e_g(2,2),e_h(2,2),e_k(2,2) 
    !将K_e划分成9个2*2的小矩阵
    area=1/2.0*(12.5**2)
    !B矩阵初始化
    B=0
    B(1,3)=12.5;B(2,6)=12.5;B(3,4)=12.5;B(3,5)=12.5
    B(1,1)=-B(1,3);B(2,2)=B(1,1);B(3,1)=B(1,1);B(3,2)=B(1,1)
    forall (i=1:3,j=1:6)
        B(i,j)=1/(2.0 * ( area) )* B(i,j)
    end forall
    !D矩阵初始化
    D=0
    D(1,1)=1
    D(2,2)=1    
    D(1,2)=p0               
    D(2,1)=p0
    D(3,3)=(1-p0)/2.0
    forall(i=1:3,j=1:3) D(i,j)=D(i,j)*E0 / (1.0-p0**2)
    K_e=matmul(matmul(transpose(B),D),B) * 1 * area
    !inv_K_e=matmul(matmul(transpose(invB),D),invB) * 1e-3 * 1/2.0 * 0.125**2
    !验证K_e与inv_K_E相等。
    !将K_e划分成9个2*2的 小矩阵
    e_a=K_e(1:2,1:2);e_b=K_e(1:2,3:4)
    e_c=K_e(1:2,5:6);e_d=K_e(3:4,1:2)
    e_e=K_e(3:4,3:4);e_f=K_e(3:4,5:6)
    e_g=K_e(5:6,1:2);e_h=K_e(5:6,3:4)
    e_k=K_e(5:6,5:6)

    !单元编号，网格划分，即确定每个单元的i,j,m。
    forall (i=1:128) E_N(i)=i
    !网格初始化
    E_i(1)=1;E_j(1)=6;E_m(1)=2
    E_i(2)=7;E_j(2)=2;E_m(2)=6
    i=3 !因为已对1、2单元赋初值，所以从3、4单元开始
    do while(i<=128)
        E_i(i)=E_i(i-2)+1
        E_j(i)=E_j(i-2)+1
        E_m(i)=E_m(i-2)+1
        if(0 == mod(E_i(i),5) ) then
            E_i(i)=E_i(i)+1
            E_j(i)=E_j(i)+1
            E_m(i)=E_m(i)+1 
        end if
        i=i+2
    end do
    i=4 !因为已对1、2单元赋初值，所以从3、4单元开始
    do while(i<=128)
        E_i(i)=E_i(i-2)+1
        E_j(i)=E_j(i-2)+1
        E_m(i)=E_m(i-2)+1
        if( mod(E_N(i)+-2,8) == 0) then
             E_i(i)=E_i(i)+1
             E_j(i)=E_j(i)+1
             E_m(i)=E_m(i)+1
        end if
        i=i+2
    end do
    open(unit=10, file='E_N_i_j_m.txt')! 打开E_N_i_j_m.txt文件,file指定文件名称。
    do i=1,128
        Element_Num(i,1)=E_i(i)
        Element_Num(i,2)=E_j(i)
        Element_Num(i,3)=E_m(i)
        write(10, "(1x,4i5)")E_N(i),E_i(i),E_j(i),E_m(i)
    end do
    El_N=transpose(Element_Num)
    !El_N(3,128)
!总刚集成
    K = 0
    forall(i=1:2)G_Element(i,i)=1
    do i = 1 , 128
        G = 0
        G_Num = 1
        do j = 1 , 3
            do a = 1 , 85
                if( El_N(j,i) == a)then
                    G((2*G_Num-1):(2*G_Num),(2*a-1):2*a ) = G_Element
                    G_Num = G_Num + 1
                    exit
                end if
            end do
        end do
        K=K + matmul(matmul(transpose(G),K_e),G)
!        do t=1,170
!            write(*, "(6i3)")(G(x,t),x=1,6)
!        end do

    end do

!   加载荷
    p(170)=-6.25
    forall (i=20:160:10) p(i)=-12.5
    open(unit=30, file='P.txt')   ! 打开P.txt文件。
    do i=1,170
        write(30,"(1x,i3,f8.2)")i,P(i)
    end do
    close(30)
!   加约束
    K(1:10,1:170)=0
    K(1:170,1:10)=0
    do i=1,10
        K(i,i)=1
        p(i)=0
    end do
    open(unit=20, file='K.txt')   ! 打开K.txt文件。
    do i=1,170
        write(20, "(170e12.3)")(K(i,j),j=1,170)
    end do
    close(20)

!求解
    data wy/170*0/
    do i=1,169
        do j=i+1,170
            temp=k(j,i)/k(i,i)
            k(j,i:170)= k(j,i:170)- k(i,i:170)*temp
            p(j)=p(j)-p(i)*temp
        end do
    end do
     do i=170,2,-1
        do j=i-1,1,-1
            temp=k(j,i)/k(i,i)
            k(j,1:170)= k(j,1:170)- k(i,1:170)*temp
            p(j)=p(j)-p(i)*temp
        end do
    end do
    forall(i=1:170)
        wy(i)=p(i)/k(i,i)
    end forall
    open(unit=40, file='W_y.txt')   ! 打开W_y.txt文件。
    write(40,"(1e12.3,/)")(wy(j),j=1,170) !输出位移数据
    close(40)

    !将wy(170)转换到N_w(85,3)
    do i=1,85
        N_w(1,i)=wy(2*i-1)
        N_w(2,i)=wy(2*i)
    end do
!    Total_N_w   Element_Num(128,3)
    El_N=transpose(Element_Num)
!    do i=1,3
!        print "(128i4)",(El_N(i,j),j=1,128)
!    end do
    do i=1,3
        do j=1,128
            do a=1,85
                if (El_N(i,j) == a)then
                    Total_N_w((2*i-1):2*i,j)=N_w(:,a)
                    cycle
                end if
            end do
        end do
    end do
    open(unit=50, file='Total_N_w.txt') ! 打开Element_Num.txt文件。
    do i=1,6
        write(50,"(128E12.3)")(Total_N_w(i,j),j=1,128)
    end do
    close(50)

    !计算应变
!    K_e=matmul(matmul(transpose(B),D),B) * 1 * area

        strain=matmul(B,Total_N_w)
    open(unit=60, file='strain.txt')   ! 打开strain.txt文件。
    do i=1,3
        write(60,"(128E12.3)")(strain(i,j),j=1,128)
    end do
    close(50)
    
    !compute stress
    stress=matmul(matmul(D,B),Total_N_w)
    open(unit=60, file='stress.txt')   ! 打开stress.txt文件。
    do i=1,3
        write(60,"(128E12.3)")(stress(i,j),j=1,128)
    end do
    close(50)
    max_stress=MAXVAL(stress)
    min_stress=MINVAL(stress)
    print "('max_stress = ',e10.3,5x,'min_stress = ',e10.3)",max_stress,min_stress
    print *,  'max_stress =  ',MAXLOC(stress)    
    print *,  'min_stress =  ',MinLOC(stress)    
end program Improving4First

