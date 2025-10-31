#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Mermaid flowchart를 이미지로 저장하는 스크립트
"""

import re
import requests
import base64
from pathlib import Path


def extract_mermaid_code(md_file: str) -> str:
    """Markdown 파일에서 mermaid 코드 블록 추출"""
    with open(md_file, 'r', encoding='utf-8') as f:
        content = f.read()
    
    # mermaid 코드 블록 찾기
    pattern = r'```mermaid\n(.*?)\n```'
    match = re.search(pattern, content, re.DOTALL)
    
    if match:
        return match.group(1).strip()
    else:
        raise ValueError("Mermaid 코드 블록을 찾을 수 없습니다.")


def export_mermaid_to_image(mermaid_code: str, output_file: str = "normalization_flowchart.png"):
    """
    Mermaid 코드를 이미지로 변환하여 저장
    
    방법 1: Mermaid.ink API 사용 (간단)
    방법 2: Playwright 사용 (로컬, 고품질)
    """
    # Mermaid 코드를 base64로 인코딩
    mermaid_b64 = base64.urlsafe_b64encode(mermaid_code.encode('utf-8')).decode('utf-8')
    
    # Mermaid.ink API URL
    api_url = f"https://mermaid.ink/img/{mermaid_b64}"
    
    print(f"🔄 Mermaid 다이어그램 다운로드 중...")
    print(f"   API: {api_url[:80]}...")
    
    try:
        # 이미지 다운로드
        response = requests.get(api_url, timeout=30)
        response.raise_for_status()
        
        # 이미지 저장
        with open(output_file, 'wb') as f:
            f.write(response.content)
        
        print(f"✅ 이미지 저장 완료: {output_file}")
        print(f"   파일 크기: {len(response.content) / 1024:.1f} KB")
        
        return output_file
        
    except Exception as e:
        print(f"❌ 에러 발생: {e}")
        print("\n💡 대안 방법:")
        print("1. https://mermaid.live 접속")
        print("2. 아래 코드를 붙여넣기")
        print("3. 다이어그램이 렌더링되면 PNG/SVG로 다운로드")
        print("\nMermaid 코드:")
        print("-" * 60)
        print(mermaid_code)
        print("-" * 60)
        return None


def main():
    """메인 함수"""
    md_file = "normalization_flowchart.md"
    output_file = "normalization_flowchart.png"
    
    print(f"📄 파일 읽는 중: {md_file}")
    
    try:
        # Mermaid 코드 추출
        mermaid_code = extract_mermaid_code(md_file)
        print(f"✅ Mermaid 코드 추출 완료 ({len(mermaid_code)} 글자)")
        
        # 이미지로 변환
        result = export_mermaid_to_image(mermaid_code, output_file)
        
        if result:
            print(f"\n🎉 완료! 이미지 파일: {Path(result).absolute()}")
        
    except FileNotFoundError:
        print(f"❌ 파일을 찾을 수 없습니다: {md_file}")
    except ValueError as e:
        print(f"❌ {e}")
    except Exception as e:
        print(f"❌ 예상치 못한 에러: {e}")


if __name__ == "__main__":
    main()

