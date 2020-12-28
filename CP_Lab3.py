from collections import namedtuple
import hashlib
import hmac
import random
import bitarray

Point = namedtuple("Point", "x y") # создание именованного списка с именем Point (x,y)
point_O = 'Оригинал'

#Parameters schema digital signature
print("***Параметры схемы ЦП:***")
p = 57896044618658097711785492504343953926634992332820282019728792003956564821041  # 256bits- module elliptic curve (E)
a = 7                    
b = 43308876546767276905765904595650931995942111794451039583252968842033849580414 #a and b are coefficients (E)
m = 57896044618658097711785492504343953927082934583725450622380973592137631069619 #order of the group of points EC (E)
q = 57896044618658097711785492504343953927082934583725450622380973592137631069619 #the order of the cyclic subgroup of the group of points of the elliptic curve E
P = Point(2, 4018974056539037503335449422937059775635739389905545080690979365213431566280)#Point P !=0
print("Модуль эллиптической кривой p -> ",p)
print("Эллиптическая кривая E, задаваемая коэффициентами a,b ∈ F_p ->",a, b)
print("Целое число m – порядок группы точек эллиптической кривой E ->", m)
print("Простое число q – порядок циклической подгруппы группы точек эллиптической кривой E ->",q)
print("Точка P≠O  эллиптической кривой E, удовлетворяющая равенству qP=O ->", P)

def valid(P):
    if P == point_O:
        return True
    else:
        return ((P.y**2 - (P.x**3 + a * P.x + b)) % p == 0 and 0 <= P.x < p and 0 <= P.y < p)


def gcd(a, b):
    s, old_s = 0, 1
    t, old_t = 1, 0
    r, old_r = b, a

    while r != 0:
        quotient = old_r // r
        old_r, r = r, old_r - quotient * r
        old_s, s = s, old_s - quotient * s
        old_t, t = t, old_t - quotient * t

    return old_r, old_s, old_t


def inverse_of(n):
    g, x, y = gcd(n, p)
    assert (n * x + p * y) % p == g, 'Ошибка с НОД'

    if g != 1:
        # Или n равно 0, или p не является простым.
        raise ValueError(
            '{} has no multiplicative inverse '
            'modulo {}'.format(n, p))
    else:
        return x % p


def add(P, Q): #2 точки проверяются на равенство с нулевой точки (на вход функции)
    zero = Point(0, 0)
    if P == zero:
        return Q
    if Q == zero:
        return P
    if P == point_O:
        return Q
    elif Q == point_O:
        return P
    elif (Q.x == P.x) and (Q.y == -P.y):
        return point_O
    else:
        if P == Q:
            с = (3 * P.x**2 + a) * inverse_of(2 * P.y)
        else:
            с = (Q.y - P.y) * inverse_of(Q.x - P.x)
        x = (с**2 - P.x - Q.x) % p
        y = (с * (P.x - x) - P.y) % p
        return Point(x, y)


def bits(n):  # Генерирует двоичные разряды n, начиная с наименее значимого бита.
    while n:
        yield n & 1
        n >>= 1


# Возвращает результат n * x, вычисленный алгоритмом удвоения-сложения
def double_and_add(n, x):
    result = Point(0, 0)
    addend = x
    for bit in bits(n):
        if bit == 1:
            result = add(addend, result)
            if result == point_O:
                return point_O
        addend = add(addend, addend)
        if addend == point_O:
            return point_O
    return result


def SGN(M, d):
    m = int(M, 16) 
    e = m % q
    if e == 0:
        e = 1
    r, s = 0, 0
    while (r == 0) or (s == 0):
        k = random.randint(1, q - 1)
        C = double_and_add(k, P)
        r = C.x % q
        s = (r * d + k * e) % q
    signature = str(format(r, 'x') + format(s, 'x'))
    len_r = r.bit_length()
    return signature, len_r


def VERIFY(signature, M, Q, len_r):
    signature = int(signature, 16)
    signature = str(format(signature, 'b'))
    r = signature[0: len_r]
    s = signature[len_r:]
    r = int(r, 2)
    s = int(s, 2)
    m = int(M, 16)
    e = m % q
    if e == 0:
        e = 1
    v = pow(e, q - 2, q)
    z1 = (s * v) % q
    z2 = ((-r) * v) % q
    C_1 = double_and_add(z1, P)
    C_2 = double_and_add(z2, Q)
    C = Point(0, 0)
    C = add(C_1, C_2)
    r2 = C.x % q
    if r == r2:
        return True
    else:
        return False


def KDF(T):
    hash_object = hashlib.sha512(T.encode(encoding="utf-8", errors="strict"))
    hex_dig = hash_object.hexdigest()
    return hex_dig


def MAC(T, k):
    k = int(k, 2)
    k = str(format(k, 'x')).encode()
    T = T.encode()
    signing = hmac.new(k, T, hashlib.sha512)
    hex_dig = signing.hexdigest()
    return hex_dig


s_A = random.randint(1, q - 1)
S_A = double_and_add(s_A, P)
Id_A = 'This is A'
s_B = random.randint(1, q - 1)
S_B = double_and_add(s_B, P)
Id_B = 'This is B'
O_I = '011001110101111100101001001010101100101001010101011101110100101'
h2 = '110110100101001101101001100100110001111010101011001001001001011'
h3 = '101010001010100101110000011101111111100101001010101100010101010'

# цифры - номера пунктов соответсвенно алгоритму


def startB(Id_A, K_a):
    # Сторона В получает от А пару (Id_A, K_A)
    if Id_A != 'This is A':
        return ValueError('Завершение протокола')
    if (not valid(K_a)) or (double_and_add(m // q, K_a) == point_O):  # 3
        return ValueError('Завершение протокола, ошибка параметра')
    print('\nСторона В получила данные')
    k_B = random.randint(1, q - 1)		# 4
    global K_B
    K_B = double_and_add(k_B, P)
    Q_AB = double_and_add(k_B * (m // q), K_a)		# 5
    # строка Id_A  в битовую строку
    global bit_arr1
    global bit_arr2
    bit_arr1 = bitarray.bitarray()
    bit_arr1.frombytes(Id_A.encode(encoding="utf-8", errors="strict"))
    # строка Id_B  в битовую строку
    bit_arr2 = bitarray.bitarray()
    bit_arr2.frombytes(Id_B.encode(encoding="utf-8", errors="strict"))
    print('Q_AB.x = ', format(Q_AB.x, 'x'))
    str_KDF = str(format(Q_AB.x, 'b') +
                  bit_arr1.to01() + bit_arr2.to01() + O_I)
    T_AB = KDF(str_KDF)		# 6
    T_AB = int(T_AB, 16)
    T_AB = str(format(T_AB, 'b'))
    global K_AB
    K_AB = T_AB[0: 255]		# 7
    print('K_AB = ', format(int(K_AB, 2), 'x'))
    M_AB = T_AB[256: 511]
    str_SGN = str(format(K_B.x, 'b') + format(K_a.x, 'b') + bit_arr1.to01())
    aut_B, len_r = SGN(str_SGN, s_B)		# 8
    print('\nСторона В сформировала код аутентификации')
    print('Код -> ', aut_B)
    str_MAC = str(h2 + format(K_B.x, 'b') + format(K_a.x, 'b') +
                  bit_arr2.to01() + bit_arr1.to01())
    tag_B = MAC(str_MAC, M_AB)
    print('Метка подтверждения B : ', tag_B)
    startA(Id_B, 'Cert_B', K_B, aut_B, tag_B, len_r)		# 9
    print('Сторона В отправила данные стороне А')


def startA(Id_B, cert, K_B, aut_B, tag_B, len_r):
    if cert != 'Cert_B':		# 10
        return ValueError('Невалидность сертификата')
    str_VERIFY = str(format(K_B.x, 'b') + format(K_A.x, 'b') + bit_arr1.to01())
    if not VERIFY(aut_B, str_VERIFY, S_B, len_r):		# 11
        return ValueError('Невозможность аутентификации')
    if (not valid(K_B)) or (double_and_add(m // q, K_B) == point_O):  # 12
        return ValueError('Завершение протокола, ошибка параметра')
    Q_BA = double_and_add(k_A * (m // q), K_B)		# 13
    print('Q_BA.x = ', format(Q_BA.x, 'x'))
    str_KDF = str(format(Q_BA.x, 'b') +
                  bit_arr1.to01() + bit_arr2.to01() + O_I)
    T_BA = KDF(str_KDF)		# 14
    T_BA = int(T_BA, 16)
    T_BA = str(format(T_BA, 'b'))
    K_BA = T_BA[0: 255]		# 15
    print('K_BA = ', format(int(K_BA, 2), 'x'))
    global M_BA
    M_BA = T_BA[256: 511]
    str_MAC = str(h2 + format(K_B.x, 'b') + format(K_A.x, 'b') +
                  bit_arr2.to01() + bit_arr1.to01())
    tag_B2 = MAC(str_MAC, M_BA)		# 16
    if tag_B != tag_B2:
        return ValueError('Невозможность аутентификации')
    print('Метка подтверждения В прошла проверку')
    str_MAC = str(h3 + format(K_A.x, 'b') + format(K_B.x, 'b') +
                  bit_arr1.to01() + bit_arr2.to01())
    tag_A = MAC(str_MAC, M_BA)		# 17
    print('Метка подтверждения А : ', tag_A)
    startB_2(tag_A)		# 18
    print('Сторона В отправила данные стороне А')


def startB_2(tag_A):
    str_MAC = str(h3 + format(K_A.x, 'b') + format(K_B.x, 'b') +
                  bit_arr1.to01() + bit_arr2.to01())
    tag_A2 = MAC(str_MAC, M_BA)
    if tag_A != tag_A2:
        return ValueError('Невозможность аутентификации')
    print('Метка подтверждения A прошла проверку')


k_A = random.randint(1, q - 1)  # 1
K_A = double_and_add(k_A, P)

startB(Id_A, K_A)  # 2
M_AB = ''
M_BA = ''
print('Удаление ключей M_AB и М_ВА')
K_AB = int(K_AB, 2)
K_AB = str(format(K_AB, 'x'))
print('Общий ключ сторон: ', K_AB)

str = str(input())
