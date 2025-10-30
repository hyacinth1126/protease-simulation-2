# 배포 가이드 (Deployment Guide)

이 프로젝트를 배포하는 방법을 안내합니다.

## 🚀 Streamlit Cloud 배포 (추천, 무료)

### 1단계: GitHub에 코드 업로드

```bash
# Git 저장소 초기화 (아직 안 했다면)
git init
git add .
git commit -m "Initial commit"

# GitHub 저장소 생성 후
git remote add origin https://github.com/YOUR_USERNAME/YOUR_REPO_NAME.git
git push -u origin main
```

### 2단계: Streamlit Cloud에서 배포

1. [share.streamlit.io](https://share.streamlit.io) 접속
2. "Sign in with GitHub" 클릭
3. "New app" 버튼 클릭
4. 다음 정보 입력:
   - **Repository**: `YOUR_USERNAME/YOUR_REPO_NAME`
   - **Branch**: `main`
   - **Main file path**: `app.py`
5. "Deploy!" 클릭

### 3단계: 배포 완료

배포 후 자동으로 생성되는 URL:
```
https://YOUR_APP_NAME.streamlit.app
```

이 URL을 공유하면 누구나 접속할 수 있습니다!

---

## 🌐 다른 클라우드 플랫폼 배포

### Heroku

1. `Procfile` 생성:
```
web: streamlit run app.py --server.port=$PORT --server.address=0.0.0.0
```

2. Heroku CLI로 배포:
```bash
heroku create your-app-name
git push heroku main
```

### AWS/Azure/GCP

각 플랫폼의 컨테이너 배포 방식을 사용하거나, Dockerfile을 생성하여 배포합니다.

---

## 💻 로컬 네트워크 접속

같은 Wi-Fi에 연결된 다른 기기에서 접속하려면:

1. 터미널에서 앱 실행:
```bash
streamlit run app.py --server.address 0.0.0.0
```

2. 터미널에 표시되는 "Network URL"을 다른 기기에서 접속:
```
http://YOUR_LOCAL_IP:8501
```

---

## 📦 필수 파일

배포를 위해 다음 파일들이 필요합니다:
- ✅ `requirements.txt` - Python 패키지 의존성
- ✅ `.streamlit/config.toml` - Streamlit 설정
- ✅ `app.py` - 메인 애플리케이션
- ✅ `analysis.py` - 분석 로직
- ✅ `plot.py` - 시각화 로직
- ✅ `fitc_peptide_timeseries.csv` - 샘플 데이터 (선택사항)

---

## ⚠️ 주의사항

1. **데이터 파일**: 샘플 데이터(`fitc_peptide_timeseries.csv`)도 함께 업로드해야 합니다.
2. **민감한 정보**: API 키나 비밀번호가 있다면 환경 변수로 관리하세요.
3. **리소스**: Streamlit Cloud 무료 플랜은 CPU/메모리 제한이 있습니다.

---

## 🔧 문제 해결

### 배포 후 에러 발생
- `requirements.txt`에 모든 필수 패키지가 있는지 확인
- 로컬에서 `pip install -r requirements.txt` 후 테스트
- Streamlit Cloud 로그 확인

### 접속 속도 느림
- 데이터 파일 크기 확인
- 불필요한 파일 제거
- 코드 최적화 검토

