import streamlit as st


def render_footer() -> None:
    """고정 하단 푸터 표시"""
    st.markdown(
        """
        <style>
        .app-footer {
            position: fixed;
            left: 0;
            right: 0;
            bottom: 0;
            background: rgba(250, 250, 250, 0.95);
            padding: 8px 16px;
            text-align: center;
            font-size: 0.9rem;
            color: #6c757d;
            border-top: 1px solid #e6e6e6;
            z-index: 100;
        }
        .app-footer a { color: inherit; text-decoration: underline; }
        </style>
        <div class="app-footer">
            Powered by <b><a href="https://github.com/hyacinth1126" target="_blank" rel="noopener noreferrer">hyacinth1126</a></b>, 2025.
        </div>
        """,
        unsafe_allow_html=True,
    )


