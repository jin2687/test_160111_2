// FLIP (Fluid-Implicit Particle) シミュレーションクラス
class FLIPSimulation {
    constructor(canvas) {
        this.canvas = canvas;
        this.ctx = canvas.getContext('2d');

        // キャンバスサイズ
        this.width = 800;
        this.height = 600;
        this.canvas.width = this.width;
        this.canvas.height = this.height;

        // シミュレーションパラメータ
        this.gravity = -9.8;
        this.dt = 0.016; // タイムステップ
        this.viscosity = 0.01;

        // グリッド設定
        this.gridSize = 10;
        this.gridWidth = Math.floor(this.width / this.gridSize);
        this.gridHeight = Math.floor(this.height / this.gridSize);

        // パーティクル
        this.particles = [];
        this.particleCount = 2000;

        // グリッド速度フィールド
        this.u = []; // x方向速度
        this.v = []; // y方向速度
        this.uPrev = [];
        this.vPrev = [];

        // マウス操作
        this.mouseX = 0;
        this.mouseY = 0;
        this.prevMouseX = 0;
        this.prevMouseY = 0;
        this.mouseDown = false;

        this.init();
        this.setupEventListeners();
    }

    init() {
        // グリッド初期化
        this.initGrids();

        // パーティクル初期化
        this.particles = [];
        const startX = this.width * 0.3;
        const startY = this.height * 0.3;
        const size = Math.sqrt(this.particleCount);
        const spacing = 3;

        for (let i = 0; i < this.particleCount; i++) {
            const row = Math.floor(i / size);
            const col = i % size;
            this.particles.push({
                x: startX + col * spacing,
                y: startY + row * spacing,
                vx: 0,
                vy: 0,
                color: `hsl(${200 + Math.random() * 40}, 100%, ${50 + Math.random() * 20}%)`
            });
        }
    }

    initGrids() {
        this.u = [];
        this.v = [];
        this.uPrev = [];
        this.vPrev = [];

        for (let i = 0; i <= this.gridWidth; i++) {
            this.u[i] = [];
            this.uPrev[i] = [];
            for (let j = 0; j <= this.gridHeight; j++) {
                this.u[i][j] = 0;
                this.uPrev[i][j] = 0;
            }
        }

        for (let i = 0; i <= this.gridWidth; i++) {
            this.v[i] = [];
            this.vPrev[i] = [];
            for (let j = 0; j <= this.gridHeight; j++) {
                this.v[i][j] = 0;
                this.vPrev[i][j] = 0;
            }
        }
    }

    setupEventListeners() {
        this.canvas.addEventListener('mousedown', (e) => {
            this.mouseDown = true;
            const rect = this.canvas.getBoundingClientRect();
            this.mouseX = e.clientX - rect.left;
            this.mouseY = e.clientY - rect.top;
            this.prevMouseX = this.mouseX;
            this.prevMouseY = this.mouseY;
        });

        this.canvas.addEventListener('mousemove', (e) => {
            if (this.mouseDown) {
                const rect = this.canvas.getBoundingClientRect();
                this.prevMouseX = this.mouseX;
                this.prevMouseY = this.mouseY;
                this.mouseX = e.clientX - rect.left;
                this.mouseY = e.clientY - rect.top;
            }
        });

        this.canvas.addEventListener('mouseup', () => {
            this.mouseDown = false;
        });

        this.canvas.addEventListener('mouseleave', () => {
            this.mouseDown = false;
        });
    }

    // パーティクルからグリッドへ速度を転送
    transferToGrid() {
        // グリッドをクリア
        for (let i = 0; i <= this.gridWidth; i++) {
            for (let j = 0; j <= this.gridHeight; j++) {
                this.u[i][j] = 0;
                this.v[i][j] = 0;
            }
        }

        // 重み付きサンプリング
        const weightU = [];
        const weightV = [];
        for (let i = 0; i <= this.gridWidth; i++) {
            weightU[i] = [];
            weightV[i] = [];
            for (let j = 0; j <= this.gridHeight; j++) {
                weightU[i][j] = 0;
                weightV[i][j] = 0;
            }
        }

        // パーティクルの速度をグリッドに転送
        for (const p of this.particles) {
            const gx = p.x / this.gridSize;
            const gy = p.y / this.gridSize;

            const i = Math.floor(gx);
            const j = Math.floor(gy);

            // バイリニア補間で転送
            if (i >= 0 && i < this.gridWidth && j >= 0 && j < this.gridHeight) {
                const fx = gx - i;
                const fy = gy - j;

                const weights = [
                    (1 - fx) * (1 - fy),
                    fx * (1 - fy),
                    (1 - fx) * fy,
                    fx * fy
                ];

                const indices = [
                    [i, j], [i + 1, j],
                    [i, j + 1], [i + 1, j + 1]
                ];

                for (let k = 0; k < 4; k++) {
                    const [gi, gj] = indices[k];
                    const w = weights[k];

                    if (gi >= 0 && gi <= this.gridWidth && gj >= 0 && gj <= this.gridHeight) {
                        this.u[gi][gj] += w * p.vx;
                        weightU[gi][gj] += w;
                        this.v[gi][gj] += w * p.vy;
                        weightV[gi][gj] += w;
                    }
                }
            }
        }

        // 重みで正規化
        for (let i = 0; i <= this.gridWidth; i++) {
            for (let j = 0; j <= this.gridHeight; j++) {
                if (weightU[i][j] > 0) this.u[i][j] /= weightU[i][j];
                if (weightV[i][j] > 0) this.v[i][j] /= weightV[i][j];
            }
        }
    }

    // グリッド上で速度場を更新（圧力投影法）
    updateGrid() {
        // 前のグリッド状態を保存
        for (let i = 0; i <= this.gridWidth; i++) {
            for (let j = 0; j <= this.gridHeight; j++) {
                this.uPrev[i][j] = this.u[i][j];
                this.vPrev[i][j] = this.v[i][j];
            }
        }

        // 重力を適用
        for (let i = 0; i <= this.gridWidth; i++) {
            for (let j = 0; j <= this.gridHeight; j++) {
                this.v[i][j] += this.gravity * this.dt;
            }
        }

        // マウス操作による力
        if (this.mouseDown) {
            const dx = this.mouseX - this.prevMouseX;
            const dy = this.mouseY - this.prevMouseY;
            const force = 50;

            const gi = Math.floor(this.mouseX / this.gridSize);
            const gj = Math.floor(this.mouseY / this.gridSize);

            for (let di = -2; di <= 2; di++) {
                for (let dj = -2; dj <= 2; dj++) {
                    const i = gi + di;
                    const j = gj + dj;
                    if (i >= 0 && i <= this.gridWidth && j >= 0 && j <= this.gridHeight) {
                        const dist = Math.sqrt(di * di + dj * dj);
                        const factor = Math.max(0, 1 - dist / 3);
                        this.u[i][j] += dx * force * factor;
                        this.v[i][j] += dy * force * factor;
                    }
                }
            }
        }

        // 境界条件
        this.applyBoundaryConditions();

        // 圧力投影（簡易版）- 非圧縮性を維持
        this.projectVelocity();
    }

    applyBoundaryConditions() {
        // 壁での速度を0に
        for (let j = 0; j <= this.gridHeight; j++) {
            this.u[0][j] = 0;
            this.u[this.gridWidth][j] = 0;
        }

        for (let i = 0; i <= this.gridWidth; i++) {
            this.v[i][0] = 0;
            this.v[i][this.gridHeight] = 0;
        }
    }

    // 圧力投影法で速度場を非圧縮にする
    projectVelocity() {
        const iterations = 20;
        const divergence = [];
        const pressure = [];

        // 配列初期化
        for (let i = 0; i <= this.gridWidth; i++) {
            divergence[i] = [];
            pressure[i] = [];
            for (let j = 0; j <= this.gridHeight; j++) {
                divergence[i][j] = 0;
                pressure[i][j] = 0;
            }
        }

        // 発散を計算
        for (let i = 1; i < this.gridWidth; i++) {
            for (let j = 1; j < this.gridHeight; j++) {
                divergence[i][j] = (
                    this.u[i + 1][j] - this.u[i][j] +
                    this.v[i][j + 1] - this.v[i][j]
                ) / this.gridSize;
            }
        }

        // ヤコビ反復で圧力を解く
        for (let iter = 0; iter < iterations; iter++) {
            const pNew = [];
            for (let i = 0; i <= this.gridWidth; i++) {
                pNew[i] = [];
                for (let j = 0; j <= this.gridHeight; j++) {
                    pNew[i][j] = pressure[i][j];
                }
            }

            for (let i = 1; i < this.gridWidth; i++) {
                for (let j = 1; j < this.gridHeight; j++) {
                    pNew[i][j] = (
                        pressure[i - 1][j] + pressure[i + 1][j] +
                        pressure[i][j - 1] + pressure[i][j + 1] -
                        divergence[i][j] * this.gridSize * this.gridSize
                    ) / 4;
                }
            }

            pressure.splice(0, pressure.length, ...pNew);
        }

        // 圧力勾配を速度から引く
        for (let i = 1; i < this.gridWidth; i++) {
            for (let j = 1; j < this.gridHeight; j++) {
                this.u[i][j] -= (pressure[i][j] - pressure[i - 1][j]) / this.gridSize;
                this.v[i][j] -= (pressure[i][j] - pressure[i][j - 1]) / this.gridSize;
            }
        }
    }

    // グリッドからパーティクルへ速度を転送（FLIP法）
    transferToParticles() {
        for (const p of this.particles) {
            const gx = p.x / this.gridSize;
            const gy = p.y / this.gridSize;

            const i = Math.floor(gx);
            const j = Math.floor(gy);

            if (i >= 0 && i < this.gridWidth && j >= 0 && j < this.gridHeight) {
                const fx = gx - i;
                const fy = gy - j;

                // バイリニア補間
                const newVX =
                    this.u[i][j] * (1 - fx) * (1 - fy) +
                    this.u[i + 1][j] * fx * (1 - fy) +
                    this.u[i][j + 1] * (1 - fx) * fy +
                    this.u[i + 1][j + 1] * fx * fy;

                const newVY =
                    this.v[i][j] * (1 - fx) * (1 - fy) +
                    this.v[i + 1][j] * fx * (1 - fy) +
                    this.v[i][j + 1] * (1 - fx) * fy +
                    this.v[i + 1][j + 1] * fx * fy;

                const oldVX =
                    this.uPrev[i][j] * (1 - fx) * (1 - fy) +
                    this.uPrev[i + 1][j] * fx * (1 - fy) +
                    this.uPrev[i][j + 1] * (1 - fx) * fy +
                    this.uPrev[i + 1][j + 1] * fx * fy;

                const oldVY =
                    this.vPrev[i][j] * (1 - fx) * (1 - fy) +
                    this.vPrev[i + 1][j] * fx * (1 - fy) +
                    this.vPrev[i][j + 1] * (1 - fx) * fy +
                    this.vPrev[i + 1][j + 1] * fx * fy;

                // FLIP: 差分を追加（0.95はFLIP、0.05はPIC）
                const flipRatio = 0.95;
                p.vx = flipRatio * (p.vx + newVX - oldVX) + (1 - flipRatio) * newVX;
                p.vy = flipRatio * (p.vy + newVY - oldVY) + (1 - flipRatio) * newVY;
            }
        }
    }

    // パーティクルを移動
    advectParticles() {
        for (const p of this.particles) {
            p.x += p.vx * this.dt * 100;
            p.y += p.vy * this.dt * 100;

            // 境界条件
            if (p.x < 0) {
                p.x = 0;
                p.vx = 0;
            }
            if (p.x > this.width) {
                p.x = this.width;
                p.vx = 0;
            }
            if (p.y < 0) {
                p.y = 0;
                p.vy = 0;
            }
            if (p.y > this.height) {
                p.y = this.height;
                p.vy = 0;
            }
        }
    }

    // シミュレーションステップ
    step() {
        this.transferToGrid();
        this.updateGrid();
        this.transferToParticles();
        this.advectParticles();
    }

    // 描画
    draw() {
        // 背景をクリア
        this.ctx.fillStyle = '#000';
        this.ctx.fillRect(0, 0, this.width, this.height);

        // パーティクルを描画
        for (const p of this.particles) {
            this.ctx.fillStyle = p.color;
            this.ctx.beginPath();
            this.ctx.arc(p.x, p.y, 2, 0, Math.PI * 2);
            this.ctx.fill();
        }
    }

    // アニメーションループ
    animate() {
        this.step();
        this.draw();
        requestAnimationFrame(() => this.animate());
    }

    reset() {
        this.init();
    }

    setGravity(value) {
        this.gravity = value;
    }

    setViscosity(value) {
        this.viscosity = value;
    }

    setParticleCount(value) {
        this.particleCount = value;
        this.init();
    }
}

// 初期化
document.addEventListener('DOMContentLoaded', () => {
    const canvas = document.getElementById('canvas');
    const simulation = new FLIPSimulation(canvas);

    // コントロール設定
    const gravityInput = document.getElementById('gravity');
    const gravityValue = document.getElementById('gravityValue');
    gravityInput.addEventListener('input', (e) => {
        const value = parseFloat(e.target.value);
        gravityValue.textContent = value.toFixed(1);
        simulation.setGravity(value);
    });

    const viscosityInput = document.getElementById('viscosity');
    const viscosityValue = document.getElementById('viscosityValue');
    viscosityInput.addEventListener('input', (e) => {
        const value = parseFloat(e.target.value);
        viscosityValue.textContent = value.toFixed(3);
        simulation.setViscosity(value);
    });

    const particleInput = document.getElementById('particleCount');
    const particleValue = document.getElementById('particleValue');
    particleInput.addEventListener('input', (e) => {
        const value = parseInt(e.target.value);
        particleValue.textContent = value;
        simulation.setParticleCount(value);
    });

    const resetButton = document.getElementById('reset');
    resetButton.addEventListener('click', () => {
        simulation.reset();
    });

    // アニメーション開始
    simulation.animate();
});
