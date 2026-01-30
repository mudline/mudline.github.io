/**************************
 * 核心工具函数
 **************************/
function linearInterpolation(xList, yList, x) {
  const n = xList.length;
  if (n === 0) throw new Error("输入列表不能为空");
  if (n !== yList.length) throw new Error("列表长度不一致");
  if (x <= xList[0]) return yList[0];
  if (x >= xList[n - 1]) return yList[n - 1];

  for (let i = 1; i < n; i++) {
    if (x <= xList[i]) {
      const xLeft = xList[i - 1];
      const xRight = xList[i];
      const yLeft = yList[i - 1];
      const yRight = yList[i];
      return yLeft + (yRight - yLeft) * (x - xLeft) / (xRight - xLeft);
    }
  }
  return yList[n - 1];
}

function newtonIteration(func, initialGuess, tol = 1e-10, maxIter = 200) {
  let x = initialGuess;
  // 数值微分（贴近Python fsolve逻辑）
  for (let i = 0; i < maxIter; i++) {
    const f = func(x);
    if (Math.abs(f) < tol) return x;
    const h = 1e-8; // 数值微分步长
    const fDeriv = (func(x + h) - func(x)) / h;
    if (Math.abs(fDeriv) < 1e-12) break; // 避免除以零
    x = x - f / fDeriv;
  }
  return x;
}

// 第一类一阶贝塞尔函数 J1(x)
function besselJ1(x) {
  const absX = Math.abs(x);
  if (absX < 1e-6) return x / 2 - x ** 3 / 16 + x ** 5 / 384;
  return Math.sin(x) / x - Math.cos(x) - (Math.sin(x) + x * Math.cos(x)) / (x ** 2) / 3;
}

// 第二类一阶贝塞尔函数 Y1(x)
function besselY1(x) {
  if (x <= 0) return -Infinity;
  const absX = Math.abs(x);
  if (absX < 1e-6) return -2 / (Math.PI * x) + x / 2 - x ** 3 / 16;
  return (-Math.cos(x) / x - Math.sin(x)) + (Math.cos(x) - x * Math.sin(x)) / (x ** 2) / 3;
}

// 第一类一阶贝塞尔函数导数 J1'(x)
function besselJ1Prime(x) {
  return (besselJ1(x - 1e-8) - besselJ1(x)) / (-1e-8); // 数值微分
}

// 第二类一阶贝塞尔函数导数 Y1'(x)
function besselY1Prime(x) {
  if (x <= 0) return -Infinity;
  return (besselY1(x - 1e-8) - besselY1(x)) / (-1e-8); // 数值微分
}

/**************************
 * 水文要素计算
 **************************/
// 线性波的弥散关系（周期 => 静水中波长）
function waveLength(T, d, g = 9.81) {
  // 定义弥散方程 f(L) = L - gT²tanh(2πd/L)/(2π)
  const func = (L) => {
    if (L === 0) return 1e10;
    return L - g * Math.pow(T, 2) * Math.tanh(2 * Math.PI * d / L) / (2 * Math.PI);
  };

  // 初始猜测与求解
  const initialGuess = 100.0;
  const LSolution = newtonIteration(func, initialGuess);

  // 浅水/深水判断
  if (d >= LSolution * 0.5) {
    return g * Math.pow(T, 2) / (2 * Math.PI);
  } else {
    return LSolution;
  }
}

// 波浪参数在水流作用下的变形计算
function waveOfCurrent(Ls, Hs, Uc, d, g = 9.81) {
  if (Uc === 0) return [Ls, Hs, Math.sqrt(2 * Math.PI * Ls / (g * Math.tanh(2 * Math.PI * d / Ls)))];

  // 定义水流中波长的求解方程
  const func = (L) => {
    if (L === 0) return 1e10;
    const k = 2 * Math.PI / L;
    const ks = 2 * Math.PI / Ls;
    const C = Uc + Math.sqrt(g * L * Math.tanh(k * d) / (2 * Math.PI));
    return L / Ls - Math.tanh(k * d) / (((1 - Uc / C) ** 2) * Math.tanh(ks * d));
  };

  // 求解水流中的波长L
  const L = newtonIteration(func, Ls);

  // 求解水流中的波高H
  const k = 2 * Math.PI / L;
  const ks = 2 * Math.PI / Ls;
  const A = 1 + 2 * k * d / Math.sinh(2 * k * d);
  const As = 1 + 2 * ks * d / Math.sinh(2 * ks * d);
  const C = Uc + Math.sqrt(g * L * Math.tanh(k * d) / (2 * Math.PI));
  const H = Math.sqrt(1 - Uc / C) * Math.sqrt(Ls / L) * Math.sqrt(As / A) * Math.pow((1 + Uc * (2 - A) / (C * A)), -0.5) * Hs;

  // 求解水流中的周期T
  const T = Math.sqrt(2 * Math.PI * L / (g * Math.tanh(k * d)));

  return [L, H, T];
}

/**************************
 * JTS 145-2015 方法
 **************************/
function loadLinearWaveJTS(D, z1, z2, H, L_wave, d, Cd, Cm, rho = 1025.0, g = 9.81) {
  const L = L_wave;
  const z_top = Math.max(z1, z2);
  const z_bot = Math.min(z1, z2);
  const gamma = rho * g; // 水的重度
  const A = 0.25 * Math.PI * (D ** 2); // 桩的横截面积

  // 规范图10.3.2-8 计算alpha
  const funcAlpha = (H, d, L) => {
    const x = d / L;

    const alpha01 = (x) => {
      let x_trunc = x;
      if (x_trunc > 0.4) x_trunc = 0.4;
      if (x_trunc < 0.1) x_trunc = 0.1;
      return 1.7 - 8.86 * x_trunc + 47.0 * (x_trunc ** 2) - 115.13 * (x_trunc ** 3) + 106.79 * (x_trunc ** 4);
    };

    const alpha02 = (x) => {
      let x_trunc = x;
      if (x_trunc > 0.35) x_trunc = 0.35;
      if (x_trunc < 0.06) x_trunc = 0.06;
      return 2.39 - 18.04 * x_trunc + 100.24 * (x_trunc ** 2) - 258.85 * (x_trunc ** 3) + 252.24 * (x_trunc ** 4);
    };

    const alpha03 = (x) => {
      let x_trunc = x;
      if (x_trunc > 0.35) x_trunc = 0.35;
      if (x_trunc < 0.03) x_trunc = 0.03;
      return 3.11 - 30.98 * x_trunc + 200.07 * (x_trunc ** 2) - 588.07 * (x_trunc ** 3) + 634.29 * (x_trunc ** 4);
    };

    const alpha04 = (x) => {
      let x_trunc = x;
      if (x_trunc > 0.25) x_trunc = 0.25;
      if (x_trunc < 0.02) x_trunc = 0.02;
      return 3.15 - 28.268 * x_trunc + 204.29 * (x_trunc ** 2) - 770.46 * (x_trunc ** 3) + 1122.54 * (x_trunc ** 4);
    };

    const alpha05 = (x) => {
      let x_trunc = x;
      if (x_trunc > 0.25) x_trunc = 0.25;
      if (x_trunc < 0.015) x_trunc = 0.015;
      return 3.08656 - 20.67846 * x_trunc + 125.184 * (x_trunc ** 2) - 451.23435 * (x_trunc ** 3) + 659.48673 * (x_trunc ** 4);
    };

    const alpha06 = (x) => {
      let x_trunc = x;
      if (x_trunc > 0.24) x_trunc = 0.24;
      if (x_trunc < 0.018) x_trunc = 0.018;
      return 3.068 - 15.64157 * x_trunc + 94.16392 * (x_trunc ** 2) - 422.79434 * (x_trunc ** 3) + 758.53721 * (x_trunc ** 4);
    };

    const alpha07 = (x) => {
      let x_trunc = x;
      if (x_trunc > 0.21) x_trunc = 0.21;
      if (x_trunc < 0.02) x_trunc = 0.02;
      return 3.02246 - 8.29662 * x_trunc + 42.57 * (x_trunc ** 2) - 349.50818 * (x_trunc ** 3) + 902.94414 * (x_trunc ** 4);
    };

    const HdRatio = H / d;
    let alpha;
    if (HdRatio <= 0.1) {
      alpha = alpha01(x);
    } else if (HdRatio > 0.1 && HdRatio <= 0.2) {
      alpha = linearInterpolation([0.1, 0.2], [alpha01(x), alpha02(x)], HdRatio);
    } else if (HdRatio > 0.2 && HdRatio <= 0.3) {
      alpha = linearInterpolation([0.2, 0.3], [alpha02(x), alpha03(x)], HdRatio);
    } else if (HdRatio > 0.3 && HdRatio <= 0.4) {
      alpha = linearInterpolation([0.3, 0.4], [alpha03(x), alpha04(x)], HdRatio);
    } else if (HdRatio > 0.4 && HdRatio <= 0.5) {
      alpha = linearInterpolation([0.4, 0.5], [alpha04(x), alpha05(x)], HdRatio);
    } else if (HdRatio > 0.5 && HdRatio <= 0.6) {
      alpha = linearInterpolation([0.5, 0.6], [alpha05(x), alpha06(x)], HdRatio);
    } else if (HdRatio > 0.6 && HdRatio <= 0.7) {
      alpha = linearInterpolation([0.6, 0.7], [alpha06(x), alpha07(x)], HdRatio);
    } else {
      alpha = alpha07(x);
    }
    return alpha;
  };

  // 规范图10.3.2-9 计算beta
  const funcBeta = (H, d, L) => {
    const x = d / L;

    const beta01 = (x) => {
      let x_trunc = x;
      if (x_trunc > 0.4) x_trunc = 0.4;
      if (x_trunc < 0.1) x_trunc = 0.1;
      return 1.52 - 4.18295 * x_trunc + 12.404 * (x_trunc ** 2) - 15.71045 * (x_trunc ** 3) + 6.97574 * (x_trunc ** 4);
    };

    const beta02 = (x) => {
      let x_trunc = x;
      if (x_trunc > 0.35) x_trunc = 0.35;
      if (x_trunc < 0.06) x_trunc = 0.06;
      return 3.77267 - 49.24525 * x_trunc + 382.20851 * (x_trunc ** 2) - 1411.98075 * (x_trunc ** 3) + 2008.29674 * (x_trunc ** 4);
    };

    const beta03 = (x) => {
      let x_trunc = x;
      if (x_trunc > 0.35) x_trunc = 0.35;
      if (x_trunc < 0.03) x_trunc = 0.03;
      return 3.61361 - 31.181 * x_trunc + 155.76364 * (x_trunc ** 2) - 365.64597 * (x_trunc ** 3) + 330.09188 * (x_trunc ** 4);
    };

    const beta04 = (x) => {
      let x_trunc = x;
      if (x_trunc > 0.25) x_trunc = 0.25;
      if (x_trunc < 0.02) x_trunc = 0.02;
      return 3.54476 - 18.15887 * x_trunc + 19.7748 * (x_trunc ** 2) + 131.02512 * (x_trunc ** 3) - 279.51738 * (x_trunc ** 4);
    };

    const beta05 = (x) => {
      let x_trunc = x;
      if (x_trunc > 0.25) x_trunc = 0.25;
      if (x_trunc < 0.015) x_trunc = 0.015;
      return 3.53367 - 5.6562 * x_trunc - 131.80141 * (x_trunc ** 2) + 787.76957 * (x_trunc ** 3) - 1247.16593 * (x_trunc ** 4);
    };

    const beta06 = (x) => {
      let x_trunc = x;
      if (x_trunc > 0.24) x_trunc = 0.24;
      if (x_trunc < 0.018) x_trunc = 0.018;
      return 3.72176 + 3.94242 * x_trunc - 280.34457 * (x_trunc ** 2) + 1568.39851 * (x_trunc ** 3) - 2636.2 * (x_trunc ** 4);
    };

    const beta07 = (x) => {
      let x_trunc = x;
      if (x_trunc > 0.21) x_trunc = 0.21;
      if (x_trunc < 0.02) x_trunc = 0.02;
      return 3.82541 + 13.86295 * x_trunc - 360.4737 * (x_trunc ** 2) + 1609.45613 * (x_trunc ** 3) - 2036.50631 * (x_trunc ** 4);
    };

    const HdRatio = H / d;
    let beta;
    if (HdRatio <= 0.1) {
      beta = beta01(x);
    } else if (HdRatio > 0.1 && HdRatio <= 0.2) {
      beta = linearInterpolation([0.1, 0.2], [beta01(x), beta02(x)], HdRatio);
    } else if (HdRatio > 0.2 && HdRatio <= 0.3) {
      beta = linearInterpolation([0.2, 0.3], [beta02(x), beta03(x)], HdRatio);
    } else if (HdRatio > 0.3 && HdRatio <= 0.4) {
      beta = linearInterpolation([0.3, 0.4], [beta03(x), beta04(x)], HdRatio);
    } else if (HdRatio > 0.4 && HdRatio <= 0.5) {
      beta = linearInterpolation([0.4, 0.5], [beta04(x), beta05(x)], HdRatio);
    } else if (HdRatio > 0.5 && HdRatio <= 0.6) {
      beta = linearInterpolation([0.5, 0.6], [beta05(x), beta06(x)], HdRatio);
    } else if (HdRatio > 0.6 && HdRatio <= 0.7) {
      beta = linearInterpolation([0.6, 0.7], [beta06(x), beta07(x)], HdRatio);
    } else {
      beta = beta07(x);
    }
    return beta;
  };

  // 规范图10.3.2-10 计算gamma_P和gamma_M
  const funcGammaPM = (d, L) => {
    const x = d / L;
    let A = 4.218;
    let B = -57.763;
    let C = 357.662;
    let D = -768.832;
    const gamma_P = A + B * x + C * x ** 2 + D * x ** 3;
    A = 4.671;
    B = -62.267;
    C = 370.421;
    D = -765.935;
    const gamma_M = A + B * x + C * x ** 2 + D * x ** 3;
    return [gamma_P, gamma_M];
  };

  // 求解分量最大力（K1-K4）
  const temp_1 = 4 * Math.PI * z_top / L;
  const temp_2 = 4 * Math.PI * z_bot / L;
  const temp_3 = Math.sinh(temp_1);
  const temp_4 = Math.sinh(temp_2);
  const temp_5 = 8 * Math.sinh(4 * Math.PI * d / L);
  const K1 = (temp_1 - temp_2 + temp_3 - temp_4) / temp_5;

  const temp1_K2 = Math.sinh(2 * Math.PI * z_top / L);
  const temp2_K2 = Math.sinh(2 * Math.PI * z_bot / L);
  const temp3_K2 = Math.cosh(2 * Math.PI * d / L);
  const K2 = (temp1_K2 - temp2_K2) / temp3_K2;

  const temp1_K3 = 1 / Math.sinh(4 * Math.PI * d / L);
  const temp2_K3 = (Math.PI ** 2 * (z_top - z_bot) ** 2) / (4 * L ** 2);
  const temp3_K3 = Math.PI * (z_top - z_bot) / (8 * L) * Math.sinh(4 * Math.PI * z_top / L);
  const temp4_K3 = 1 / 32 * (Math.cosh(4 * Math.PI * z_top / L) - Math.cosh(4 * Math.PI * z_bot / L));
  const K3 = temp1_K3 * (temp2_K3 + temp3_K3 - temp4_K3);

  const temp1_K4 = 1 / Math.cosh(2 * Math.PI * d / L);
  const temp2_K4 = 2 * Math.PI * (z_top - z_bot) / L * Math.sinh(2 * Math.PI * z_top / L);
  const temp3_K4 = Math.cosh(2 * Math.PI * z_top / L) - Math.cosh(2 * Math.PI * z_bot / L);
  const K4 = temp1_K4 * (temp2_K4 - temp3_K4);

  // 最大速度分力和最大惯性分力
  let P_D_max = Cd * gamma * D * (H ** 2) * K1 / 2.0;
  let P_I_max = Cm * gamma * A * H * K2 / 2.0;

  // z_bot底部的力矩
  let M_D_max = Cd * gamma * D * (H ** 2) * L * K3 / (2 * Math.PI);
  let M_I_max = Cm * gamma * A * H * L * K4 / (4 * Math.PI);

  // 计算修正系数
  const alpha = funcAlpha(H, d, L);
  const beta = funcBeta(H, d, L);

  // 系数修正逻辑
  if ((H / d <= 0.2 && d / L >= 0.2) || (H / d > 0.2 && d / L >= 0.35)) {
    // 无需修正
  } else if ((H / d <= 0.2 && d / L < 0.2) || (H / d > 0.2 && d / L < 0.35)) {
    let gamma_P = 1.0, gamma_M = 1.0;
    if (0.04 <= (d / L) && (d / L) <= 0.2) {
      [gamma_P, gamma_M] = funcGammaPM(d, L);
    }
    P_D_max = P_D_max * alpha;
    P_I_max = P_I_max * gamma_P;
    M_D_max = M_D_max * beta;
    M_I_max = M_I_max * gamma_M;
  }

  // 小尺度/大尺度桩判断，合成最终结果
  let P_max, M_max;
  if (D / L <= 0.2) {
    if (P_D_max <= (0.5 * P_I_max)) {
      P_max = P_I_max;
      M_max = M_I_max;
    } else {
      P_max = P_D_max * (1 + 0.25 * Math.pow(P_I_max / P_D_max, 2));
      M_max = M_D_max * (1 + 0.25 * Math.pow(M_I_max / M_D_max, 2));
    }
  } else {
    const x = D / L;
    const CM = 1.98924 +
                2.20959 * x +
                (-3.93435) * (x ** 2) +
                (-92.09772) * (x ** 3) +
                308.76789 * (x ** 4) +
                (-365.1082) * (x ** 5) +
                150.64515 * (x ** 6);
    P_I_max = CM * gamma * A * H * K2 / 2.0;
    M_I_max = CM * gamma * A * H * L * K4 / (4 * Math.PI);
    P_max = P_I_max;
    M_max = M_I_max;
  }

  return [P_max, M_max];
}

/**************************
 * NBT 11084-2023 方法
 **************************/
// 海流力计算
function loadLinearCurrentNBT(D, z1, z2, u, Cd, rho = 1028.0) {
  const f_current = 0.5 * rho * Math.pow(u, 2);
  const F_current = Cd * f_current * D * Math.abs(z2 - z1);
  const M_current = F_current * 0.5 * Math.abs(z2 - z1);
  return [F_current, M_current];
}

// 波浪力计算
function loadLinearWaveBook(D, z1, z2, H, L_wave, d, Cd, Cm, rho = 1028.0, g = 9.81) {
  [z1, z2] = [Math.min(z1, z2), Math.max(z1, z2)];
  const gamma = rho * g; // 水的重度

  // 小尺度桩（D/L <= 0.2）
  if (D / L_wave <= 0.2) {
    const k = 2 * Math.PI / L_wave;

    // 计算K1-K4
    const K1 = (2 * k * (z2 - z1) + Math.sinh(2 * k * z2) - Math.sinh(2 * k * z1)) / (8 * Math.sinh(2 * k * d));
    const K2 = (Math.sinh(k * z2) - Math.sinh(k * z1)) / Math.cosh(k * d);
    const K3 = 1 / (32 * Math.sinh(2 * k * d)) * (
      2 * Math.pow(k, 2) * Math.pow(z2 - z1, 2) +
      2 * k * (z2 - z1) * Math.sinh(2 * k * z2) -
      (Math.cosh(2 * k * z2) - Math.cosh(2 * k * z1))
    );
    const K4 = 1 / Math.cosh(k * d) * (
      k * (z2 - z1) * Math.sinh(k * z2) -
      (Math.cosh(k * z2) - Math.cosh(k * z1))
    );

    // 最大速度分力与惯性分力
    const P_D_max = Cd * gamma * D * Math.pow(H, 2) * K1 / 2.0;
    const P_I_max = Cm * gamma * Math.PI * Math.pow(D, 2) * H * K2 / 8.0;

    // 底部力矩
    const M_D_max = Cd * gamma * D * Math.pow(H, 2) * L_wave * K3 / (2.0 * Math.PI);
    const M_I_max = Cm * gamma * Math.pow(D, 2) * H * L_wave * K4 / 16.0;

    // 合成结果
    let P_max, M_max;
    if (P_D_max <= (0.5 * P_I_max)) {
      P_max = P_I_max;
      M_max = M_I_max;
    } else {
      P_max = P_D_max * (1 + 0.25 * Math.pow(P_I_max / P_D_max, 2));
      M_max = M_D_max * (1 + 0.25 * Math.pow(M_I_max / M_D_max, 2));
    }
    return [P_max, M_max];
  }

  // 大尺度桩（D/L > 0.2）
  else {
    const k = 2.0 * Math.PI / L_wave;
    const ka = k * D * 0.5;

    // 计算Aka
    const J1p = besselJ1Prime(ka);
    const Y1p = besselY1Prime(ka);
    const Aka = 1.0 / Math.sqrt(Math.pow(J1p, 2) + Math.pow(Y1p, 2));

    // 最大水平力
    const P_max = 2.0 * rho * g * H / (k * Math.cosh(k * d)) * Aka * (1 / k) * (
      Math.sinh(k * z2) - Math.sinh(k * z1)
    );

    // 最大力矩（简化公式）
    const M_max = P_max * Math.abs(z2 - z1) * 0.5;

    return [P_max, M_max];
  }
}