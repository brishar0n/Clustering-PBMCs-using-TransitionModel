import streamlit as st
import numpy as np
import pandas as pd
import plotly.express as px
import time

# Set Dashboard Configuration
st.set_page_config(
    page_title="Advanced PBMC Simulation",
    layout="wide",
    page_icon="🩺",
)

# Sidebar for Navigation
st.sidebar.title("Navigation")
menu = st.sidebar.radio(
    "Go to",
    ["Home", "Simulation"]
)

# Extended Cluster Data
clusters = ["T-cells", "B-cells", "Monocytes",
            "Macrophages", "Dendritic Cells"]
default_transition_matrix = np.array(
    [[0.6, 0.2, 0.1, 0.05, 0.05],  # T-cells
     [0.1, 0.7, 0.1, 0.05, 0.05],  # B-cells
     [0.15, 0.15, 0.5, 0.1, 0.1],  # Monocytes
     [0.05, 0.05, 0.3, 0.5, 0.1],  # Macrophages
     [0.1, 0.1, 0.2, 0.2, 0.4]]    # Dendritic Cells
)

# Initialize Session State
if "transition_matrix" not in st.session_state:
    st.session_state.transition_matrix = default_transition_matrix.copy()
if "distributions" not in st.session_state:
    st.session_state.distributions = None
if "current_distribution" not in st.session_state:
    st.session_state.current_distribution = np.array(
        [1.0, 0.0, 0.0, 0.0, 0.0])  # Start with all T-cells
if "simulation_running" not in st.session_state:
    st.session_state.simulation_running = False

# Distinct colors for clusters
cluster_colors = {
    "T-cells": "blue",
    "B-cells": "green",
    "Monocytes": "orange",
    "Macrophages": "purple",
    "Dendritic Cells": "red"
}

# ----------------------------------------------
# Home Page
# ----------------------------------------------
if menu == "Home":
    st.title("Advanced PBMC Simulation Dashboard")
    st.markdown("""
    ### About the Project
    The **Advanced PBMC Simulation Dashboard** is a tool for simulating the dynamic behavior of 
    peripheral blood mononuclear cells (PBMCs) over time. PBMCs include several immune cell types, 
    such as T-cells, B-cells, Monocytes, Macrophages, and Dendritic Cells, which are critical 
    for the immune response.

    ### Features
    - **Real-Time Simulation**: Observe dynamic transitions between extended cell types.
    - **Interactive Controls**: Adjust transition probabilities, initial distributions, and speed.
    - **Advanced Visualizations**: Analyze results using line charts, bar charts, and pie charts.
    """)

# ----------------------------------------------
# Simulation Page
# ----------------------------------------------
elif menu == "Simulation":
    st.title("PBMC Simulation")
    st.markdown(
        "Configure the simulation parameters below and observe real-time transitions between PBMC cell types.")

    # Input Controls
    st.sidebar.subheader("Simulation Controls")
    time_steps = st.sidebar.slider("Number of Time Steps", 1, 100, 30)
    speed = st.sidebar.slider("Speed (seconds per step)", 0.1, 2.0, 0.5)
    stochasticity = st.sidebar.slider("Stochastic Factor (%)", 0, 50, 10) / 100
    initial_distribution = st.sidebar.radio(
        "Select Initial Distribution of Cells:",
        ["All in T-cells", "All in B-cells",
            "All in Monocytes", "Custom Distribution"]
    )

    # Set initial distributions
    if initial_distribution == "All in T-cells":
        st.session_state.current_distribution = np.array(
            [1.0, 0.0, 0.0, 0.0, 0.0])
    elif initial_distribution == "All in B-cells":
        st.session_state.current_distribution = np.array(
            [0.0, 1.0, 0.0, 0.0, 0.0])
    elif initial_distribution == "All in Monocytes":
        st.session_state.current_distribution = np.array(
            [0.0, 0.0, 1.0, 0.0, 0.0])
    else:
        t_cells = st.sidebar.slider("T-cells (%)", 0.0, 1.0, 0.2)
        b_cells = st.sidebar.slider("B-cells (%)", 0.0, 1.0, 0.2)
        monocytes = st.sidebar.slider("Monocytes (%)", 0.0, 1.0, 0.2)
        macrophages = st.sidebar.slider("Macrophages (%)", 0.0, 1.0, 0.2)
        dendritic_cells = 1.0 - (t_cells + b_cells + monocytes + macrophages)
        st.session_state.current_distribution = np.array(
            [t_cells, b_cells, monocytes, macrophages, dendritic_cells]
        )

    st.write("Initial Distribution of Cells:",
             st.session_state.current_distribution)

    # Buttons to Start/Stop Simulation
    col1, col2 = st.columns(2)
    with col1:
        if st.button("Run Simulation"):
            st.session_state.simulation_running = True
    with col2:
        if st.button("Stop Simulation"):
            st.session_state.simulation_running = False

    # Real-Time Simulation
    if st.session_state.simulation_running:
        st.markdown("### Simulation Running...")

        # Placeholders for Visualizations
        line_chart_placeholder = st.empty()

        # Initialize Distributions List
        distributions = [st.session_state.current_distribution]

        for step in range(time_steps):
            if not st.session_state.simulation_running:
                break

            # Apply Transition Matrix with Stochasticity
            transition_matrix = st.session_state.transition_matrix + np.random.normal(
                0, stochasticity, st.session_state.transition_matrix.shape
            )
            # Ensure no negative probabilities
            transition_matrix = np.maximum(transition_matrix, 0)
            transition_matrix = transition_matrix / \
                transition_matrix.sum(axis=1, keepdims=True)

            # Compute Next Distribution
            next_distribution = np.dot(distributions[-1], transition_matrix)
            distributions.append(next_distribution)

            # Update Line Chart
            df = pd.DataFrame(distributions, columns=clusters)
            df["Time Step"] = range(len(distributions))
            line_chart_fig = px.line(
                df, x="Time Step", y=clusters,
                title="Cluster Proportions Over Time (Line Chart)",
                color_discrete_map=cluster_colors,
                labels={"value": "Proportion", "variable": "Cluster"}
            )
            line_chart_placeholder.plotly_chart(
                line_chart_fig, use_container_width=True)

            # Pause for Next Step
            time.sleep(speed)

        # Save Distributions
        st.session_state.distributions = distributions
        st.success("Simulation Completed!")

        # Show Results in a Modal
        with st.expander("Click to view detailed results"):
            st.write("### Final Distribution")
            final_distribution = pd.DataFrame({
                "Cluster": clusters,
                "Proportion": distributions[-1]
            })
            st.dataframe(final_distribution)

            # Pie Chart
            st.write("### Final Distribution (Pie Chart)")
            pie_chart_fig = px.pie(
                final_distribution, names="Cluster", values="Proportion",
                title="Final Cluster Distribution",
                color="Cluster",
                color_discrete_map=cluster_colors
            )
            st.plotly_chart(pie_chart_fig)

    elif not st.session_state.simulation_running:
        st.warning("Simulation Stopped. Click 'Run Simulation' to restart.")
