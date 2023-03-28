# 분석 환경 모듈
import numpy as np
import pandas as pd
import psycopg2 as pg
import time
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import style

# 데이터 분석 모듈
from CoolProp.CoolProp import PropsSI
from sklearn.preprocessing import MinMaxScaler, QuantileTransformer, PolynomialFeatures
from sklearn.metrics import mean_absolute_error, mean_squared_error
from sklearn.linear_model import LinearRegression
from scipy.integrate import ode

# 내가 만든 모듈 임포트
from Model import *
from PhysicalProperty import *
from StructuredQuery import *
from Numeric import *

style.use('seaborn-talk')
krfont = {'family':'Times New Roman', 'weight':'bold', 'size':10}
matplotlib.rc('font', **krfont)
matplotlib.rcParams['axes.unicode_minus'] = False

if __name__ == "__main__":
    # 클래스 정의 및  인스턴스 생성
    pro = PhysicalProperty()
    mod = ModelOSV()
    sql = StructuredQuery()

    # CHF 데이터베이스 연결
    load_chf_query = "SELECT * FROM rawdata_chf_tb"
    chf_tb = sql.read_sql(load_chf_query, db_engine)

    # Rawdata에서 Physical property 추가 (for water)
    chf_tb['pcrit'] = round(chf_tb[['refri']].apply(lambda x: PropsSI(x[0],'pcrit') * 10**-5, axis=1),6) # [bar]
    chf_tb['tsat'] = round(chf_tb[['p', 'refri']].apply(lambda x: PropsSI('T', 'P', x[0] * 1e5, 'Q', 0, x[1]), axis=1),6) # [K]
    chf_tb['kf'] = round(chf_tb[['p', 'refri']].apply(lambda x: PropsSI('L', 'P', x[0] * 1e5, 'Q', 0, x[1]), axis=1),6) # [W/mK]
    chf_tb['kv'] = round(chf_tb[['p', 'refri']].apply(lambda x: PropsSI('L', 'P', x[0] * 1e5, 'Q', 1, x[1]), axis=1),6) # [W/mK]
    chf_tb['muf'] = round(chf_tb[['p', 'refri']].apply(lambda x: PropsSI('V', 'P', x[0] * 1e5, 'Q', 1, x[1]), axis=1),12) #  [Pa s]
    chf_tb['muv'] = round(chf_tb[['p', 'refri']].apply(lambda x: PropsSI('V', 'P', x[0] * 1e5, 'Q', 1, x[1]), axis=1),12) # [Pa s]
    chf_tb['hfo'] = round(chf_tb[['p', 'refri']].apply(lambda x: PropsSI('H', 'P', x[0] * 1e5, 'Q', 0, x[1]), axis=1),6) # [J/kgK]
    chf_tb['hgo'] = round(chf_tb[['p', 'refri']].apply(lambda x: PropsSI('H', 'P', x[0] * 1e5, 'Q', 1, x[1]), axis=1),6) #[J/kgK]
    chf_tb['lam'] = round(chf_tb['hgo'] - chf_tb['hfo'],6) # [J/kgK]
    chf_tb['rhof'] = round(chf_tb[['p', 'refri']].apply(lambda x: PropsSI('D', 'P', x[0] * 1e5, 'Q', 0, x[1]), axis=1),6) #[kg/m3]
    chf_tb['rhov'] = round(chf_tb[['p', 'refri']].apply(lambda x: PropsSI('D', 'P', x[0] * 1e5, 'Q', 1, x[1]), axis=1),6) # [kg/m3]
    chf_tb['v'] = round(chf_tb['g'] / chf_tb['rhof'],6) # [m/s]
    chf_tb['cpf'] = round(chf_tb[['p', 'refri']].apply(lambda x: PropsSI('C', 'P', x[0] * 1e5, 'Q', 0, x[1]), axis=1),6) # [J/kgK]
    chf_tb['cpv'] = round(chf_tb[['p', 'refri']].apply(lambda x: PropsSI('C', 'P', x[0] * 1e5, 'Q', 1, x[1]), axis=1),6) # [J/kgK]
    chf_tb['sigma'] = round(chf_tb[['p', 'refri']].apply(lambda x: PropsSI('I', 'P', x[0] * 1e5, 'Q', 0, x[1]), axis=1),6) # [N/m]
    chf_tb['ti'] = round(chf_tb['tsat'] - chf_tb['dtin'] - 273.15,6) # [K]

    # Rawdata에서 Dimensionless number 추가
    chf_tb['xi'] = round(chf_tb[['cpf', 'dtin', 'lam']].apply(lambda x: POSV.Xi(x[0], x[1], x[2]), axis = 1),6)
    chf_tb['xout'] = round(chf_tb[['q', 'doi', 'dio', 'geo', 'hsur', 'g', 'cpf', 'lam', 'dtin', 'lh']].apply(lambda x: POSV.Xout(x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9]), axis=1),6)
    chf_tb['de'] = round(chf_tb[['doi', 'dio', 'geo', 'hsur', 'dh']].apply(lambda x: POSV.De(x[0], x[1], x[2], x[3], x[4]), axis=1),6)
    chf_tb['pe'] = round(chf_tb[['dh', 'g', 'cpf', 'kf']].apply(lambda x: POSV.Pe(x[0], x[1], x[2], x[3]), axis=1),6)
    chf_tb['re'] = round(chf_tb[['g', 'dh', 'muf']].apply(lambda x: POSV.Re(x[0], x[1], x[2]), axis=1),6)
    chf_tb['we'] = round(chf_tb[['rhof', 'v', 'dh', 'sigma']].apply(lambda x: POSV.We(x[0], x[1], x[2], x[3]), axis=1),6)
    chf_tb['bd'] = round(chf_tb[['rhof', 'rhov', 'dh', 'sigma']].apply(lambda x: POSV.Bd(x[0], x[1], x[2], x[3]), axis=1),6)
    chf_tb['bo'] = round(chf_tb[['q', 'lam', 'g']].apply(lambda x: POSV.Bo(x[0], x[1], x[2]), axis = 1),6)
    chf_tb['bo_el'] = round(chf_tb[['q', 'rhof', 'rhov', 'sigma', 'v', 'lam']].apply(lambda x: POSV.Bo_el (x[0], x[1], x[2], x[3], x[4], x[5]), axis=1),6)
    chf_tb['ga'] = round(chf_tb[['rhof', 'dh', 'muf']].apply(lambda x: POSV.Ga(x[0], x[1], x[2]), axis=1),6)
    chf_tb['gz'] = round(chf_tb[['g', 'cpf', 'kf', 'dh']].apply(lambda x: POSV.Gz(x[0], x[1], x[2], x[3]), axis=1),6)
    chf_tb['z'] = round(chf_tb[['muf', 'rhof', 'dh', 'sigma']].apply(lambda x: POSV.Z(x[0], x[1], x[2], x[3]), axis=1),6)
    chf_tb['fr'] = round(chf_tb[['v', 'dh']].apply(lambda x: POSV.Fr(x[0], x[1]), axis=1),6)
    chf_tb['ca'] = round(chf_tb[['muf', 'v', 'sigma', 'rhof']].apply(lambda x: POSV.Ca(x[0], x[1], x[2], x[3]), axis=1),6)
    chf_tb['pr'] = round(chf_tb[['cpf', 'muf', 'kf']].apply(lambda x: POSV.Pr(x[0], x[1], x[2]), axis=1),6)

    print("chf_tb에 Physical property 계산 완료")

    # OSV, CHF 설정 모델 삽입 문구
    print('chf correlation을 선택하시오: ')
    print('(1) Saha and Zuber  (2) Modified Saha and Zuber  (3) Unal')
    print('(4) Levy                         (5) Bowring                                    (6) Jeong')
    m_idx = int(input('-----너의 선택은? : '))

    # Xt 초기값 계산하는 알고리즘
    """
    res_tb: chf_tb (rawdata table)에서 뽑아낸 결과만 저장하는 테이블

    """
    # 결과테이블 새로 생성
    res_tb = chf_tb

    # def sub_xt_return(p, pcrit, dh, g, lam, m_idx):
    if m_idx == 1:
        for i, row in res_tb.iterrows():
            res_tb.loc[i, 'dt_sz'], res_tb.loc[i, 'x_sz'] = MOSV.SahaZuber(res_tb.loc[i, 'q'], res_tb.loc[i, 'rhof'], res_tb.loc[i, 'dh'], res_tb.loc[i, 'g'],
                                                res_tb.loc[i, 'cpf'], res_tb.loc[i, 'kf'], res_tb.loc[i, 'pe'], res_tb.loc[i, 'lam'])
            res_tb.loc[i,'st_cal'] = res_tb.loc[i,'q']* 10 ** 6 / (res_tb.loc[i,'g'] * res_tb.loc[i,'cpf'] * res_tb.loc[i,'dt_sz'])
        print("Saha and Zuber correlation에 대한 St, Xosv 계산이 모두 끝났습니다.")
    elif m_idx == 2:
        for i, row in res_tb.iterrows():
            res_tb.loc[i, 'dt_psz'], res_tb.loc[i, 'x_psz'] = MOSV.ParkSahaZuber(res_tb.loc[i, 'q'], res_tb.loc[i, 'rhof'], res_tb.loc[i, 'dh'], res_tb.loc[i, 'g'],
                                                res_tb.loc[i, 'cpf'], res_tb.loc[i, 'kf'], res_tb.loc[i, 'pe'], res_tb.loc[i, 'lam'])
            res_tb.loc[i,'st_cal'] = res_tb.loc[i,'q'] * 10 ** 6 / (res_tb.loc[i,'g'] * res_tb.loc[i,'cpf'] * res_tb.loc[i,'dt_psz'])
        print("Park, Saha and Zuber correlation에 대한 St, Xosv 계산이 모두 끝났습니다.")
    elif m_idx == 3:
        for i, row in res_tb.iterrows():
            res_tb.loc[i, 'dt_unal'], res_tb.loc[i, 'x_unal'] = MOSV.Unal(res_tb.loc[i, 'q'], res_tb.loc[i, 'pr'], res_tb.loc[i, 'dh'], res_tb.loc[i, 'v'],
                                                res_tb.loc[i, 'cpf'], res_tb.loc[i, 'kf'], res_tb.loc[i, 're'], res_tb.loc[i, 'refri'], res_tb.loc[i, 'lam'])
            res_tb.loc[i,'st_cal'] = res_tb.loc[i,'q'] * 10 ** 6 / (res_tb.loc[i,'g'] * res_tb.loc[i,'cpf'] * res_tb.loc[i,'dt_unal'])
        print("Unal correlation에 대한 St, Xosv 계산이 모두 끝났습니다.")
    elif m_idx == 4:
        for i, row in res_tb.iterrows():
            res_tb.loc[i, 'dt_levy'], res_tb.loc[i, 'x_levy'] = MOSV.Levy(res_tb.loc[i, 'sigma'], res_tb.loc[i, 'dh'], res_tb.loc[i, 'rhof'], res_tb.loc[i, 'muf'],
                                                res_tb.loc[i, 'kf'], res_tb.loc[i, 're'], res_tb.loc[i, 'pr'], res_tb.loc[i, 'cpf'],
                                                res_tb.loc[i, 'g'], res_tb.loc[i, 'q'], res_tb.loc[i, 'lam'], res_tb.loc[i, 'v'])
            res_tb.loc[i,'st_cal'] = res_tb.loc[i,'q']* 10 ** 6 / (res_tb.loc[i,'g'] * res_tb.loc[i,'cpf'] * res_tb.loc[i,'dt_levy'])
        print("Levy Model에 대한 St, Xosv 계산이 모두 끝났습니다.")
    elif m_idx == 5:
        for i, row in res_tb.iterrows():
            res_tb.loc[i, 'dt_bowr'], res_tb.loc[i, 'x_bowr'] = MOSV.Bowring(res_tb.loc[i, 'p'], res_tb.loc[i, 'q'], res_tb.loc[i, 'v'], res_tb.loc[i, 'lam'],
                                                    res_tb.loc[i, 'cpf'])
            res_tb.loc[i,'st_cal'] = res_tb.loc[i,'q'] * 10 ** 6 / (res_tb.loc[i,'g'] * res_tb.loc[i,'cpf'] * res_tb.loc[i,'dt_bowr'])
        print("Bowring correlation에 대한 St, Xosv 계산이 모두 끝났습니다.")
    else:
        for i, row in res_tb.iterrows():
            res_tb.loc[i, 'dt_js'], res_tb.loc[i, 'x_js'] = MOSV.Jeong(res_tb.loc[i, 'q'], res_tb.loc[i, 'rhof'], res_tb.loc[i, 'dh'], res_tb.loc[i, 'v'],
                                            res_tb.loc[i, 'cpf'], res_tb.loc[i, 'kf'], res_tb.loc[i, 'pe'], res_tb.loc[i, 'lam'], res_tb.loc[i,'pr'], res_tb.loc[i, 'we'])
            res_tb.loc[i,'st_cal'] = res_tb.loc[i,'q'] * 10 ** 6 / (res_tb.loc[i,'g'] * res_tb.loc[i,'cpf'] * res_tb.loc[i,'dt_js'])
        print("Jeong and Shim correlation에 대한 St, Xosv 계산이 모두 끝났습니다.")    

        # Boundary conditoin설정
    """
    - 국부조건 가설을 사용한 CHF를 계산하는 과정
    dXt       | h_b - h_l,sat   |
    --- =  |----------- |*Xe
    dXe      |         \lambda |
    """
    # 시간기록 시작 (나중에 여기에 로그도 남겨야함)
    startTime = time.time()

    # Rate eqn. 설정
    def dif(X, Xt):
        global nparam_xosv
        dXt=1+(Xt-X)/(nparam_xosv*(1-Xt));
        return dXt

    # Validation set와  error table 만드는 작업 진행
    val_chf_tb = res_tb[['source','q','xi','xe','xout']]
    val_chf_err_tb = res_tb[['source','q', 'xi','xe','xout']]

    try:
        for i, row in val_chf_tb.iterrows():
            # alpha와 gamma 계산
            # 나중에 CHF 모델 선택도 가능하게 추가
            val_chf_tb.loc[i,'alpha'] = 1.669 - 6.544 * ((res_tb.loc[i,'p']/res_tb.loc[i,'pcrit'])-0.448)**2
            val_chf_tb.loc[i,'gamma'] = 0.06523 + (0.1045 / np.sqrt(2* np.pi*(np.log(res_tb.loc[i,'p']/res_tb.loc[i,'pcrit']))**2))* np.exp(-5.413*((np.log(res_tb.loc[i,'p']/res_tb.loc[i,'pcrit'])+0.4537)**2/(np.log(res_tb.loc[i,'p']/res_tb.loc[i,'pcrit']))**2))
            
            # 앞서 계산한 OSV 상관식으로 계산한 st이 삽입된 xosv_cal 계산
            val_chf_tb.loc[i,'xosv_cal'] = -(res_tb.loc[i,'q'] * 10**6)/ (res_tb.loc[i,'st_cal'] * res_tb.loc[i,'g'] * res_tb.loc[i,'lam'])
            
            # Xt 초기값 계산 (ODE_45)
            nparam_xosv = val_chf_tb.loc[i,'xosv_cal']
            xsol, ysol = rk_ode45(dif,res_tb.loc[i, 'xi'],[val_chf_tb.loc[i,'xosv_cal'],res_tb.loc[i,'xe']],0.0001)
            val_chf_tb.loc[i,'xt_ini'] = ysol[-1]
            
            if val_chf_tb.loc[i, 'xt_ini'] < 0:
                val_chf_tb.loc[i,'xt_cal'] = 0
            else:
                val_chf_tb.loc[i,'xt_cal'] = ysol[-1]
            
            cnt = 1
            while 1:
                """
                differential rate eqn.을 돌려보니, xe = xt
                """
                # CHF heat flux를 계산 (Deng 식)
                val_chf_tb.loc[i,'term_deng'] = np.sqrt(res_tb.loc[i,'g']*val_chf_tb.loc[i,'xt_cal']*(1+val_chf_tb.loc[i,'xt_cal']**2)**3)
                val_chf_tb.loc[i,'q_cal'] = (val_chf_tb.loc[i,'alpha']/np.sqrt(res_tb.loc[i,'dh']))*np.exp(-val_chf_tb.loc[i,'gamma']*val_chf_tb.loc[i,'term_deng'])
                print('heat flux of calculation is {}'.format(val_chf_tb.loc[i, 'q_cal']))
                
                # 오류값 계산 (나중에 바꿔줄 필요 있음)
                val_chf_tb.loc[i, 'q_err'] = (val_chf_tb.loc[i,'q_cal'] - res_tb.loc[i,'q']) / res_tb.loc[i,'q']
                
                # Error
                # print('error is {}'.format(err))
                
                # 조건 판별 (break 또는 증감)
                try:
                    cnt = cnt + 1               
                    print('{}번째 시도에서 xt_cal은 {}입니다. '.format(cnt, val_chf_tb.loc[i,'xt_cal']))
                    
                    if cnt >= 100:
                        raise Exception('divergence error')
                    
                    if  np.abs(val_chf_tb.loc[i, 'q_err']) < 0.05:
                        print('{}th iteration is converged., xt_cal is {} and error of heat flux is {}.'.format(i, val_chf_tb.loc[i,'xt_cal'], val_chf_tb.loc[i,'q_err']))
                        val_chf_err_tb.loc[i,'idx_err'] = 'Converged'
                        break
                    else:
                        if val_chf_tb.loc[i, 'q_err'] > 0:
                            val_chf_tb.loc[i,'xt_cal'] = val_chf_tb.loc[i,'xt_cal'] + 0.01
                            continue
                        #elif val_chf_tb.loc[i,'xt_cal'] > 1:
                        #    val_chf_tb.loc[i,'xt_cal'] = val_chf_tb.loc[i,'xt_cal'] / 2
                        #    continue
                        else:
                            val_chf_tb.loc[i,'xt_cal'] = val_chf_tb.loc[i,'xt_cal'] - 0.01
                            continue
                    
                except Exception as e:
                    print('{}th iteration is diverged.'.format(i))
                    val_chf_err_tb.loc[i,'idx_err'] = 'Diverged'
                    break
                                                                                                        
    except Exception as e:
        print('The index {} is error occured.'.format(e))
        print(e)
    finally:
        print("{} th process is succees!".format(i))

    endTime = time.time() - startTime
    print('전체 시뮬레이션 시간은 {} 모델에 대해 데이터 개수는 {}, 총 걸린 시간은 {} ms입니다.'.format(m_idx, len(chf_tb), endTime))

    # 데이터를 PostgreSQL로 저장하기
    sql.write_sql(val_chf_tb, 'val_res_js_tb',db_engine)
    sql.write_sql(val_chf_err_tb, 'val_res_js_err_tb', db_engine)