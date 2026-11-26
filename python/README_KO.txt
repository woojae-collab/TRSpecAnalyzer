TRSpecAnalyzer - Python판 1.1.0
===============================

MATLAB 프로그램 TRSpecAnalyzer의 Python판입니다. 시간분해 분광 데이터
(TA, TR-PL, TR-Raman, TR-IR)를 분석하며 MATLAB 없이 실행됩니다.
1.1.0은 MATLAB판 v1.1.0의 창을 모두 갖추고 있고 계산 방식도 같습니다.
버전 번호는 "기능이 같은 MATLAB판의 번호"입니다. 창의 사용법은 MATLAB판
사용 설명서를 그대로 보시면 됩니다.

   1. 들어 있는 것
   2. 설치 (처음 한 번): Windows, macOS
   3. 처음 해 볼 것 (데모 데이터)
   4. 업데이트
   5. 그림에서의 마우스 동작
   6. MATLAB판과 다른 점
   7. MATLAB판과 숫자가 다를 수 있는 곳
   8. 확인한 방법
   9. 문제가 생기면 / 지우기 / 다른 운영체제 / 폴더 구성
  10. 변경 기록


1. 들어 있는 것
---------------
- 메인 창: 데이터 불러오기, 2D 맵 / 스펙트럼 / kinetics, 십자선(클릭, 드래그),
  Scale (Linear / Log / Split), 표시 범위, 컬러맵, Z 범위, 단위, Pin
- 보정: Background Corr, Chirp Corr, Crop, Set t=0, Mask, Reset all
- 여러 파일 불러와 평균 (Load & Average, Auto-align 포함)
- Export Corrected Data, Save / Load Session
- 데이터 도구(둘째 줄 버튼): Subtract data (용매나 기판만 측정한 데이터 빼기,
  scale과 time shift는 직접 입력 또는 추정), Anisotropy (불러온 parallel 데이터
  옆에 perpendicular 데이터를 불러와 parallel / perpendicular / magic angle /
  anisotropy r 중 골라 보기), ΔA <-> ΔT/T 변환. 파일의 숫자가 무엇인지
  (ΔA, ΔT/T, Intensity)는 불러올 때 고릅니다.
- 분석: Fit trace (한 파장, 또는 대역의 적분이나 평균), SVD, Global Analysis
  (coherent artifact 항, stretched, IRF FWHM 처리까지 MATLAB판과 같음),
  Target Analysis, LDA, MCR-ALS, Coherence (Time-frequency 포함, 용매만 측정한
  데이터로 용매의 vibrational coherence 빼기 포함), Auto Analysis
- 분석 창의 결과는 창을 닫아도 남고 세션 파일에 저장됩니다. 메인 창의 데이터가
  바뀌면 지워집니다.


2. 설치 (처음 한 번): Windows, macOS
------------------------------------
[Windows]
(1) Python 설치
    https://www.python.org/downloads/ 에서 Python 3.13 (64-bit)을 받습니다.
    설치 첫 화면에서 "Add python.exe to PATH"에 체크한 뒤 "Install Now".
    (3.12도 됩니다. 막 나온 최신 버전은 필요한 패키지가 아직 없을 수 있어
     3.13을 권합니다.)

(2) 이 폴더의 install_windows.bat 을 더블클릭
    이 폴더 안에 .venv 라는 폴더를 만들고 필요한 패키지 네 가지
    (numpy, scipy, matplotlib, PySide6)를 받습니다.
    인터넷 연결이 필요하고 몇 분 걸립니다. 디스크는 약 500 MB를 씁니다.
    PC의 다른 곳은 건드리지 않습니다.

(3) run_windows.bat 을 더블클릭하면 프로그램이 뜹니다.
    함께 뜨는 검은 창은 프로그램의 메시지를 보여 주는 창입니다.
    프로그램을 닫으면 같이 닫힙니다.

"Windows의 PC 보호" 창이 뜨면 "추가 정보" -> "실행"을 누르면 됩니다
(.bat 파일을 처음 실행할 때 나올 수 있습니다).

[macOS]
(1) https://www.python.org/downloads/ 에서 macOS용 Python 3.13을 받아 설치합니다.

(2) 이 폴더의 install_mac.command 를 더블클릭
    터미널 창이 열리고, 이 폴더 안에 .venv 폴더를 만들어 패키지 네 가지를
    받습니다 (인터넷 필요, 몇 분, 약 500 MB).

(3) run_mac.command 를 더블클릭하면 프로그램이 뜹니다.
    함께 열리는 터미널 창은 프로그램의 메시지를 보여 줍니다. 프로그램을 닫은 뒤에
    닫으면 됩니다.

처음에는 macOS가 인터넷에서 받은 .command 파일을 열지 않겠다고 할 수 있습니다.
그럴 때는 둘 중 하나로 하세요.
  - 시스템 설정 > 개인정보 보호 및 보안에서 아래로 내려가, 그 파일에 대한
    "그래도 열기"를 누릅니다 (macOS 14 이하에서는 파일을 오른쪽 클릭 > "열기"도
    됩니다).
  - 터미널에서: "cd " (뒤에 한 칸)를 입력하고 이 폴더를 터미널 창으로 끌어다
    놓은 뒤 Return, 그다음
        bash install_mac.command
    실행할 때는
        bash run_mac.command
"적절한 접근 권한이 없어서 실행할 수 없다"고 나오면 터미널에서 이 폴더로 이동해
        chmod +x install_mac.command run_mac.command

macOS용 두 파일은 Linux에서만 돌려 보았습니다 (같은 셸을 씁니다). Mac에서는
실행해 보지 못했습니다.


3. 처음 해 볼 것 (데모 데이터)
------------------------------
demo 폴더의 README_demo.txt 에 있는 순서를 그대로 해 볼 수 있습니다.
(GitHub 저장소를 통째로 받은 경우 demo 폴더는 python 폴더 안이 아니라 옆에
있습니다.) 앞의 아홉 단계를 줄여 적으면:

 1. Load...            demo\demo_TA.csv -> Load
 2. Background Corr... 평균할 점 개수 10 -> Apply & Close
 3. Chirp Corr...      Auto-place -> Fit & Apply
 4. Crop...            t range -1.5 to 2800 -> Apply Crop
 5. Mask λ...          From 494, To 506 -> Add region
 6. SVD...             특이값 세 개가 나머지보다 훨씬 큽니다
 7. Global Analysis... Components 3, 시작값 1 / 10 / 100,
                       IRF t_0 와 IRF FWHM 의 "Fixed" 해제 -> Run Fit
                       기대값: tau = 0.60, 15.0, 701 ps, FWHM 0.063 ps,
                       RMS 0.060
 8. Target Analysis... Model "Sequential, all decay:  A -> B -> C -> GS",
                       "τ from GA" -> Run fit
                       기대값: 7번과 같은 tau (0.60, 15.0, 701 ps), RMS 0.060
                       (7번을 먼저 해야 합니다. t_0 와 FWHM을 거기서 가져옵니다.)
 9. Coherence...       Method "Use GA residual", t_min 0.15, t_max 4
                       -> Run Coherence
                       기대값: 150, 250 cm-1 두 피크. 150 cm-1 은 556, 606 nm
                       근처에서 최대, 580 nm 에서 마디 (7번의 residual을 씁니다)

Auto Analysis (2-9번을 한 번에): demo\demo_TA.csv 를 새로 불러온 상태에서
   Auto Analysis... -> Run
   기대값 (Report 탭): Background 처음 43개 delay, Mask 없음 (500 ... 502 nm는
   "noisier than the neighbors, but kept"), Chirp -0.27 ps (450 nm) ~ +0.56 ps
   (750 nm), Crop t = -1.7 ... 2730 ps, artifact 항 not needed (3~4 %),
   성분 3개 tau = 0.604, 15.0, 700 ps, IRF FWHM 0.063 ps, RMS 0.069,
   진동 150, 249 cm-1 (마디 580, 654 nm 근처).
   "Apply corrections to the data" -> "Send model to Global Analysis" -> 거기서
   "Run Fit" 하면 같은 값이 나옵니다. (2-5번의 손 보정과는 보정이 달라서 7번과
   숫자가 조금 다릅니다: background 43개, 다른 crop, mask 없음.)

여러 파일 평균: Load... 에서 demo\scans 의 파일 다섯 개를 한꺼번에 선택합니다.

Target Analysis 전용 데모 두 가지 (README_demo.txt 의 뒷부분에 순서가 있습니다)
- demo\demo_target.csv       갈라지는 모델 A -> B, A -> C. "equal:" 두 줄을
                             넣으면 k1..k4 = 2.99, 6.02, 40.0, 799 ps, RMS 0.052
- demo\demo_equilibrium.csv  양방향 단계 A <-> B. "zero: B 700 750" 한 줄을
                             넣으면 kf, kb, k2, k3 = 4.99, 15.0, 100.0, 1998 ps,
                             RMS 0.053

1.1.0에서 추가된 기능용 데모 두 쌍 (README_demo.txt 의 "Solvent demo", "Anisotropy demo"
부분에 순서가 있습니다)
- demo\demo_solution.csv     cyclohexane에 녹인 용질과 cyclohexane만 측정한 데이터.
  demo\demo_solvent.csv      Coherence 창의 "Solvent coherence...": A = 0.700,
                             time shift 0.0030 ps. 용매의 801 cm-1 선이 없어지고
                             215, 590 cm-1가 남습니다. 그 뒤 메인 창의 "View"에서
                             용매가 빠진 데이터를 볼 수 있습니다.
- demo\demo_solution_raw.csv   같은 두 파일의 측정된 그대로의 상태 (background,
  demo\demo_solvent_raw.csv    chirp). solution 파일을 보정한 뒤 "Solvent
                             coherence..."에서 solvent 파일을 불러오면 같은 보정이
                             적용됩니다. A = 0.699, time shift 0.0030 ps.
- demo\demo_aniso_par.csv    한 시료를 두 편광으로 측정한 데이터. "Anisotropy...",
  demo\demo_aniso_perp.csv   tail matching: G = 1.15. 이른 delay에서 r = 0.38과
                             -0.12이고 48 ps로 사라집니다.
- demo\demo_aniso_raw_par.csv   같은 두 파일의 측정된 그대로의 상태 (background, pump
  demo\demo_aniso_raw_perp.csv  산란, chirp 포함). parallel 파일을 전처리한 뒤
                             perpendicular 파일을 불러오면 같은 보정이 적용됩니다.
- demo\demo_aniso_raw_perp_late.csv  perpendicular 측정을 time zero가 0.05 ps 늦은
                             상태로 한 번 더 한 것. "Anisotropy..."의 "Estimate":
                             time shift -0.050 ps.

Time-resolved emission 데모 (Signal = Intensity로 불러옵니다. README_demo.txt 의
"Solvation demo" 부분에 순서가 있습니다)
- demo\demo_solvation.csv    용매가 완화되는 동안 띠가 장파장으로 이동하는 데이터.
                             420 ~ 760 nm 대역을 Fit trace (Integral, 성분 1개,
                             시작값 1000)로 피팅하면 형광 수명 3987 ps만 나옵니다.
                             "Band position...", Log-normal fit, 0.2 ps부터 지수
                             3개: 0.245, 2.92, 14.8 ps (만들 때 0.25, 3, 15).


4. 업데이트
-----------
새 압축을 풀어 나온 trspecgui, trspec 폴더와 README 파일들을 기존 폴더에
덮어쓰면 됩니다. .venv 는 그대로 두면 다시 설치할 필요가 없습니다. (새 폴더에
풀었다면 거기서 install_windows.bat 을 한 번 실행하면 됩니다.) 이전 판에서 저장한
세션 파일과 설정은 그대로 쓸 수 있습니다. 새 패키지가 필요한 판은 이 자리에
그렇게 적습니다.


5. 그림에서의 마우스 동작
-------------------------
- 마우스 휠: 포인터 위치를 중심으로 확대 / 축소. Shift를 누르거나 포인터를 x축
  숫자 위에 두면 x축만, Ctrl을 누르거나 y축 숫자 위에 두면 y축만 바뀝니다.
- 오른쪽 버튼으로 드래그: 그 사각형 영역으로 확대.
- 왼쪽 버튼으로 드래그: 확대된 그림을 옮깁니다. 클릭에 이미 뜻이 있는 그림
  (십자선이 있는 2D 맵, 표시선이 있는 맵과 스펙트럼)에서는 휠을 누른 채로
  드래그하세요.
- 더블클릭: 원래 범위로. 클릭에 뜻이 있는 그림에서는 오른쪽 클릭 -> "Reset view".
  휠을 끝까지 반대로 돌려도 원래 범위가 됩니다.
- 확대한 그림은 프로그램이 다시 그려도 그대로입니다: 다른 delay / 파장 선택, Pin
  추가, 같은 범위에서의 보정. 다른 데이터를 불러오거나 "Reset all"을 누르거나,
  단위 / Scale (Linear, Log, Split) / 그 축의 표시 범위를 바꾸면 원래 범위로
  돌아갑니다.
- Split 보기에서 kinetics 그림의 두 부분은 y 범위를 같이 쓰고, 맵의 두 패널은 파장
  범위를 같이 씁니다. 분석 창의 데이터 / 피팅 / residual 맵 세 개는 함께 확대됩니다.
- 범례는 끌어서 옮길 수 있고, 옮긴 자리는 그림이 다시 그려져도 유지됩니다.
- 그림 위 오른쪽 클릭: 그림 복사 / 저장, 범례 숨기기 / 자동 위치로, "Reset view",
  그리고 위 내용을 정리해 보여 주는 "Zoom with the mouse...".
- 트랙패드 (macOS): 두 손가락 스크롤이 휠, 두 손가락 클릭이 오른쪽 클릭입니다.
  "휠을 누른 채로 드래그"는 휠이 달린 마우스가 있어야 합니다.


6. MATLAB판과 다른 점
---------------------
- 세션 파일은 .trsession 입니다. MATLAB판의 .mat 세션과 서로 열 수 없습니다.
- 컬러맵 parula 대신 viridis 가 들어 있습니다.
- 확대 / 이동은 MATLAB이 그림 위에 띄우는 도구 모음 대신 마우스로 합니다 (5번).
  data tip은 없습니다. MATLAB판은 코드상 메인 창 스펙트럼을 다시 그릴 때마다 파장
  범위를 되돌리는데, 여기서는 확대가 유지됩니다.
- Split 보기의 kinetics: 오른쪽 부분에는 세로축 숫자가 없고 (두 부분이 y 범위를
  같이 씁니다), 범례는 넓은 오른쪽 부분에 놓입니다.
- Split 모드에서 맵의 위 칸에서 아래 칸으로 끌면 십자선이 따라옵니다.
- 창 모양은 항상 밝은 색입니다 (Windows의 어두운 모드를 따르지 않습니다).
  그림의 글꼴은 어느 PC에서나 같습니다 (matplotlib에 들어 있는 DejaVu Sans).
- 작은 화면이나 큰 화면 배율에서: 제목과 축 이름은 자리에 맞지 않으면 글자가 조금
  작아지거나 짧은 표현으로 바뀌고, 버튼 줄은 다음 줄로 넘어가며, 창이 낮으면 설정
  칸과 그림이 스크롤됩니다. 메인 창 오른쪽 위의 상태 글은 너무 길면 "..."으로
  줄어들고, 마우스를 올리면 전체가 보입니다.
- Target Analysis: "Export:" 글자가 버튼들 위가 아니라 같은 줄 앞에 있습니다.
  창이 낮으면 scheme 칸과 parameter 표가 조금 낮아집니다 (결과 창 자리 확보).
- Auto Analysis: Report의 긴 줄은 다음 줄로 이어집니다. 창이 좁으면 작은 그림 두
  개의 legend와 kinetics 위의 안내 글이 짧은 말로 바뀝니다. log 축 눈금은 10의
  거듭제곱 대신 보통 숫자(0.03)로 씁니다.
- 피팅 중에는 "Please wait" 창이 뜹니다 (Global Analysis는 수 초에서 수십 초).
- 설정(마지막 폴더, 단위, 컬러맵 등)은 Windows에서는
  %LOCALAPPDATA%\TRSpecAnalyzer\prefs_python.json, macOS에서는
  ~/Library/Preferences/TRSpecAnalyzer/prefs_python.json 에 저장됩니다.
- 피팅과 계산 방식은 MATLAB판과 같습니다.


7. MATLAB판과 숫자가 다를 수 있는 곳
------------------------------------
계산 방식은 같습니다. 다만 데이터가 정해 주지 못하는 값은 탐색이 지나가는 경로에
따라 달라지는데, 그 경로가 두 프로그램에서 똑같지는 않습니다.

- Auto Analysis
    * "three terms lower the RMS by ... %" (시험 삼아 넣어 보는 artifact 항):
      demo_TA에서 경로에 따라 3~4 %. "not needed"라는 판단은 같습니다.
    * artifact가 없는 데이터에 "Coherent artifact: On"을 켠 경우: 해가 가까이 두 개
      있고, 작은 차이로 어느 쪽에서 끝나는지가 갈립니다. demo_TA, 기본 Mask rule:
      양쪽 모두 IRF FWHM 0.084 ps. Mask rule "Every noisy channel": 이 판은
      0.080 ps (RMS 0.0575), MATLAB 코드는 0.084 ps (RMS 0.0569). 시간 상수
      차이는 0.6 % 이하입니다.
    * 성분이 하나 많은 피팅(양쪽 다 기각)은 기각 사유가 다를 수 있습니다.
    * rise가 없는 데이터(time zero 뒤에서 시작)인데 time zero가 있다고 본 경우:
      합성 시험 데이터 하나에서, 이 판은 탐색이 IRF 폭을 0으로 몰고 갔고 MATLAB
      코드는 그러지 않았습니다. 그러면 폭이 고정되고 (Report에 적힘) 피팅이
      나빠집니다. 그런 데이터에서 Report에 "IRF FWHM ... (fixed)"가 보이면 앞쪽
      delay를 확인하거나 잘라 주세요.
- Target Analysis: 결과 창에 "[not determined]"로 표시되는 rate constant는 다른
  숫자가 나올 수 있습니다 (README_demo.txt 의 target 데모 3번 "71 % to B",
  equilibrium 데모 3번 "4.0, 50 and 124 ps"). 데이터가 정해 주는 것(lifetime, RMS,
  제약을 넣은 뒤의 값)은 같습니다. 같은 이유로 README_demo.txt 8번에 적힌 안내
  ("The direct fit lost the IRF width ...")는 여기서는 나오지 않습니다. 탐색이 같은
  답에 바로 도달하기 때문입니다.
- Coherence의 "Per-pixel exp fit": 시간 상수를 자유롭게 두는 피팅이라, 일부
  파장에서는 다른 국소해에서 끝날 수 있습니다.
- MCR-ALS의 "Random" 시작: MATLAB과 다른 난수를 씁니다. 같은 seed를 넣어도
  MATLAB판과는 다른 시작점이 됩니다. 이 판 안에서는 같은 seed가 같은 결과를 줍니다.


8. 확인한 방법
--------------
- 창마다 MATLAB 프로그램의 코드(GNU Octave에서 실행)와 비교했습니다. 글, 표, 안내
  문구는 글자까지 같습니다. 피팅 값은 7자리 이상 (Global / Target Analysis, 두
  탐색이 같은 해에서 끝날 때), LDA와 MCR-ALS 결과는 10자리 이상, Coherence와
  Time-frequency의 맵은 12자리 이상, Auto Analysis의 피팅 표는 같은 보정
  데이터에서 6자리 이상 같습니다.
- demo\README_demo.txt 의 순서를 따라 하면 거기 적힌 숫자가 나옵니다 (3번).
- 만드는 동안 저자가 Windows에서 써 보았습니다.
- 확인하지 않은 것: macOS와 Linux의 화면 (Linux에서 검사용으로 화면 없이만
  돌렸습니다), 그리고 macOS용 두 파일.


9. 문제가 생기면 / 지우기 / 다른 운영체제 / 폴더 구성
-----------------------------------------------------
- 설치나 실행이 안 되면 검은 창에 나온 글을 복사해 두세요.
- 프로그램 안에서 오류 창이 뜨면 "Show Details"의 내용을 복사해 두세요.
- 지우기: 이 폴더를 통째로 지우면 됩니다. (설정 파일은 6번에 적은 위치에 있습니다.)
- Linux, 또는 macOS에서 직접 입력할 때:
      python3 -m venv .venv
      .venv/bin/python -m pip install -r requirements.txt
      .venv/bin/python -m trspecgui
- 폴더:  trspecgui\  화면(GUI) 코드,  trspec\  계산 코드,
         demo\  데모 데이터와 설명


10. 변경 기록
-------------
1.1.0  (2026-11-26)
   MATLAB판 1.1.0의 새 기능들을 같은 계산 방식으로 넣었습니다
   (자세한 내용은 저장소의 CHANGELOG.md).
   - 신호 종류: 파일의 숫자가 무엇인지(ΔA, ΔT/T, Intensity)를 불러올 때 고르고,
     축 이름, 컬러바, 내보내기 파일 이름이 이를 따릅니다. ΔA <-> ΔT/T 변환은
     정확한 관계식으로 합니다 ("ΔA <-> ΔT/T..." 버튼).
   - "Subtract data...": 같은 축에서 측정한 기준(용매, 기판)을 뺍니다:
     데이터 - scale x 기준(t - shift). scale과 shift는 직접 입력하거나 최소제곱으로
     추정합니다. 기준 파일은 측정한 그대로 불러오면 되고, crop과 Auto Analysis에서
     적용한 보정을 따라갑니다.
   - "Anisotropy...": 불러온 (parallel) 데이터 옆에 perpendicular 데이터를 불러와
     똑같이 보정합니다. "View"에서 parallel, perpendicular, magic angle,
     anisotropy r 중 하나를 고르면 보이는 것이 분석 대상이 됩니다. G는 직접
     입력하거나 tail matching으로 구합니다. tail matching에 제안되는 범위는
     delay의 마지막 한 decade(마지막 delay의 1/10부터 끝까지)입니다. tail
     matching은 모든 파장에 대해 G 하나를 구합니다. 창의 세 번째 그림이 tail
     범위에서 파장 하나하나의 G를 불확도와 함께 보여 주어, G 하나로 충분한지
     확인할 수 있습니다 (보기 전용이며 계산에는 쓰이지 않습니다).
     두 측정의 time zero가 다를 때는 "Time shift of the perpendicular data"에
     값을 입력하거나 "Estimate"로 구합니다. "Estimate"는 신호가 올라오는 구간에서
     perpendicular 데이터가 parallel 데이터를 가장 잘 따라가는 shift를 찾습니다
     (파장마다 배율과 offset을 따로 둡니다). 값과 함께 파장 사이의 흩어짐에서 구한
     표준편차를 알려 줍니다.
   - 분석 창의 설정과 결과가 view마다 따로 보관됩니다. view를 바꾸면 분석 창은
     닫히고, 고른 view의 결과가 다시 돌아옵니다 (세션 파일에도 저장됩니다).
     어떤 view의 데이터가 바뀌면 그 view에 보관된 것은 지워집니다: 보정, crop,
     Reset all은 전부, 새 G나 perpendicular 데이터의 새 time shift는 Parallel view를
     뺀 전부, 새 threshold는 r view만.
   - Coherence 창의 "Solvent coherence...": 용매만 측정한 데이터의 vibrational
     coherence를 FFT 전에 잔차에서 뺍니다. 배율과 time shift는 용매 모드가 있는
     대역에서 Fourier 변환의 실수부와 허수부로 추정합니다
     (M. Liebel, 학위논문, Oxford 2014).
     여기서 빼는 동안("Apply") 메인 창의 "View"에서 용매가 빠진 데이터를 고를 수
     있습니다: "Solvent coherence subtracted"(분석의 time window 안에서 용매의
     진동 부분만) 또는 "Solvent signal subtracted"(여기서 구한 배율과 shift로
     모든 delay에서 용매 신호 전체). 2D 맵, export, 모든 분석이 보이는 view의
     데이터를 씁니다. view를 바꾸면 분석 창은 닫히고, 설정과 결과는 view마다
     보관되어 그 view로 돌아오면 다시 나옵니다.
   - Coherence 창, 빼기 전과 후: 용매 빼기를 하면 빼기 전의 결과도 함께 보관됩니다.
     그림 위의 줄에서 어느 맵을 볼지 고르고, 곡선 그림에는 빼기 전의 곡선이 회색으로
     함께 그려지며, "Side by side..."는 두 맵을 같은 색 범위로 나란히 보여 주는 창을
     엽니다. export는 그림에 보이는 것을 저장합니다.
   - Coherence 창의 "Color range": |FFT|^2 맵의 색 범위(z scale)를 직접 입력할 수
     있습니다. 예를 들어 최댓값을 작게 주면 약한 모드가 잘 보입니다. "Auto"를 누르면
     자동 범위로 돌아갑니다.
   - Fit trace: 한 파장 대신 여러 파장의 대역을 피팅할 수 있습니다 ("Band, to:").
     spectral axis에 대한 적분 또는 평균을 씁니다. Anisotropy r view에서는 두 편광의
     대역 값으로 r을 계산합니다.
   - "Band position..." (새 창): 띠가 delay마다 어디에 있는지를 구합니다 (dynamic
     Stokes shift 같은 spectral shift). 위치는 first moment, 꼭대기를 지나는
     포물선의 최댓값, log-normal 띠 피팅의 최댓값 중 하나로 구하며 (창에서 선택),
     cm-1 또는 nm로 나타냅니다. nm당 intensity는 cm-1당 신호로 변환합니다. 위치의
     시간 변화는 지수 1 ~ 4개와 상수로 피팅합니다. 창에는 한 delay의 띠와 모델,
     위치와 피팅, C(t), 띠의 폭이 그려지고, 표로 내보낼 수 있습니다.
   - delay 방향의 보간: chirp 보정과 time shift("Subtract data...", perpendicular
     데이터)는 delay 간격이 일정한 구간 안에서는 cubic spline으로, 그 밖에서는 선형
     보간으로 delay 사이의 값을 구합니다. 1.0.2까지 chirp 보정은 전체가 선형이어서
     빠른 진동이 약해졌습니다(한 주기에 delay 8개: 파워가 평균 9%, 파장에 따라
     0 ~ 14%). 그래서 chirp 보정을 거친 데이터는 1.0.2와 조금 다릅니다.
     스크립트용: trspec.chirp에 delay_runs, spline_curv, interp_runs, chirp_shift.
   - 이 기능들의 세부 사항: tail matching은 두 데이터의 잡음을 모두 고려하는 직교
     회귀로 G를 구합니다. Anisotropy r view에서는 맵의 자동 색 범위가 -0.5 ... 0.5를 넘지
     않습니다. 메인 창에서 빼기를 적용하거나 바꾸거나 없애면 Coherence 창의 용매
     coherence 빼기는 꺼집니다.
   - 이 기능들을 위한 데모 데이터: demo_solution.csv와 demo_solvent.csv
     (측정된 그대로의 같은 쌍 demo_solution_raw.csv, demo_solvent_raw.csv 포함),
     demo_aniso_par.csv와 demo_aniso_perp.csv (측정된 그대로의 같은 쌍
     demo_aniso_raw_*.csv와, time zero가 0.05 ps 늦은 perpendicular 측정
     demo_aniso_raw_perp_late.csv 포함), 그리고 demo_solvation.csv
     (스펙트럼이 이동하는 time-resolved emission, 대역 피팅용). 순서와 기대되는
     숫자는 README_demo.txt 에 있습니다.
   - IRF 폭을 피팅하는 경우 (Fit trace, Global Analysis, Target Analysis):
     FWHM이 delay 범위의 절반을 넘으면, delay 간격으로 정해지지 않는 폭과
     마찬가지로 "폭을 잃었다"고 보고 두 단계(IRF 고정, 그다음 풀기)로 다시
     피팅합니다. 전에는 시작값보다 훨씬 느린 decay(ns decay에 시작값 1 ps)에서
     FWHM이 delay 범위까지 커지고 시간상수가 무의미한 값으로 끝날 수 있었습니다.
     다시 해도 그대로이면 정해지지 않은 값으로 표시됩니다.
   - Auto Analysis: Run 뒤의 상태 줄이, 정하지 못한 단계가 있으면 "No correction
     was needed" 대신 그 단계 이름을 알려 줍니다("Not determined: chirp").
   - 그림: Global Analysis 창의 map 세 개(Experimental, Global fit, Residual)
     모두에 신호 이름이 붙은 color bar가 있습니다. Global fit map에는 전에
     color bar가 없었습니다. 확대한 log y축에서 여섯 글자 이상이 되는 숫자는 공통 인수
     (x 10^n)로 적어서, 축 이름이 그림 밖으로 밀려나지 않습니다. MCR-ALS 창의
     convergence 그림은 LOF의 변화가 10배보다 작을 때 LOF 축에 숫자를 그대로
     (1.37, 1.38) 적습니다. 전에는 이 숫자들이 왼쪽 끝에서 잘렸습니다.
   - 그림의 행과 열: 그림이 여러 개인 모든 창에서, 같은 행의 그림은 위아래 끝이,
     같은 열의 그림은 왼쪽 끝이 정확히 같습니다. color bar가 없는 그림은 map의
     color bar 자리를 비워 두지 않습니다. Coherence, Target, Lifetime Density
     Analysis 창, Global Analysis 창의 map 아래 행들, Auto Analysis 창의 Residual
     탭에서 이런 그림의 오른쪽 여백은 y축 이름 왼쪽의 빈 공간과 비슷한 너비입니다.
     다음 그림은 위나 아래의 map과 너비가 정확히 같게 남습니다: Lifetime Density
     Analysis 창의 map 아래 스펙트럼(map과 파장 축이 같음), Time-Frequency 창의
     map 아래 스펙트럼, SVD 창의 그림들, Auto Analysis 창 Corrections 탭의 그림들.
     Components 탭의 DADS와 EADS는 너비가 같습니다.
   Coherence 창의 per-pixel fit 잔차는 설정이 바뀔 때까지 보관됩니다. 이전 버전의
   세션은 그대로 열립니다. 스크립트용으로 trspec.coherence에 coh_spectrum,
   solvent_vc, solvent_shifted, solvent_band, fit_solvent_coherence가 추가되었습니다.

1.0.2  (2026-11-22)
   번호를 붙인 첫 판입니다. 미리보기 12의 프로그램에, 기능이 같은 MATLAB판의
   번호를 붙였습니다. 화면과 보고서의 영어 글은 MATLAB판 1.0.2와 마찬가지로 미국식
   철자를 씁니다. macOS용 install_mac.command, run_mac.command 포함.

미리보기 판 (2026년 11월)
   12  데이터를 불러오기 전의 메인 창, Split 보기: 맵의 두 패널 사이에 있던 숫자와
       "Wavelength (nm)"를 없앰.
   11  Auto Analysis가 마지막 실행 결과를 기억: 창을 다시 열면 돌아오고
       ("[restored]"), 세션 파일에 저장됨.
   10  모든 그림에서 마우스로 확대 / 축소 / 이동.
    9  메인 창의 상태 글이 첫 줄에 남음.
    8  Auto Analysis. 고친 것: Coherence, Time-frequency, Auto-align의 등간격 delay
       격자가 마지막 delay를 반올림 오차만큼 지나서 끝날 수 있던 문제.
    7  Coherence, Time-frequency.
    6  LDA, MCR-ALS.
    5  Target Analysis.
    4  끌어 옮긴 범례에 관한 수정.
    3  범례를 끌어 옮기기, Windows 어두운 모드에서 체크박스 테두리, 작은 화면에
       맞춘 글씨.
 1, 2  메인 창, 보정, Load & Average, 세션, Fit trace, SVD, Global Analysis.


License: MIT (LICENSE 파일).  저자: 김우재, ULS Lab, 연세대학교 화학과.
코드는 Claude Code로 작성했습니다.
