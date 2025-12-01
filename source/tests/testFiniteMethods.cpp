/*
 * FiniteMethods unit tests
 */

#define private public
#define protected public
#include "../plugins/components/FiniteMethods.h"
#undef private
#undef protected

#include "../../source/gtest/gtest/gtest.h"
#include "../kernel/simulator/Simulator.h"
#include <vector>

class FiniteMethodsFixture {
public:
	FiniteMethodsFixture() : simulator(), model(simulator.getModelManager()->newModel()), finite(model) {}
	Simulator simulator;
	Model* model;
	FiniteMethods finite;
};

TEST(FiniteMethods, DiffusionNoneKeepsState) {
	FiniteMethodsFixture fx;
	fx.finite.setDiffusionModel(FiniteMethods::DiffusionModel::NONE);
	std::vector<double> current{1.0, 2.0, 3.0};
	std::vector<double> next;
	fx.finite._applyDiffusion(current, next);
	EXPECT_DOUBLE_EQ(next[0], 1.0);
	EXPECT_DOUBLE_EQ(next[1], 2.0);
	EXPECT_DOUBLE_EQ(next[2], 3.0);
}

TEST(FiniteMethods, ForwardEulerExplicitStep) {
	FiniteMethodsFixture fx;
	fx.finite.setDiffusionModel(FiniteMethods::DiffusionModel::CONSTANT);
	fx.finite.setSolver(FiniteMethods::Solver::EXPLICIT);
	fx.finite.setDiscretization(FiniteMethods::Discretization::FORWARD_EULER);
	fx.finite.setDiffusionCoefficient(1.0);
	fx.finite.setSpaceStep(1.0);
	fx.finite.setTimeStep(0.5); // alpha = 0.5
	std::vector<double> current{1.0, 0.0, 0.0};
	std::vector<double> next;
	fx.finite._applyDiffusion(current, next);
	EXPECT_DOUBLE_EQ(next.front(), 1.0);
	EXPECT_DOUBLE_EQ(next[1], 0.5);
	EXPECT_DOUBLE_EQ(next.back(), 0.0);
}

TEST(FiniteMethods, CrankNicolsonExplicitBlendsPredictor) {
	FiniteMethodsFixture fx;
	fx.finite.setDiffusionModel(FiniteMethods::DiffusionModel::CONSTANT);
	fx.finite.setSolver(FiniteMethods::Solver::EXPLICIT);
	fx.finite.setDiscretization(FiniteMethods::Discretization::CRANK_NICOLSON);
	fx.finite.setDiffusionCoefficient(1.0);
	fx.finite.setSpaceStep(1.0);
	fx.finite.setTimeStep(0.5); // alpha = 0.5
	std::vector<double> current{0.0, 1.0, 0.0};
	std::vector<double> next;
	fx.finite._applyDiffusion(current, next);
	EXPECT_DOUBLE_EQ(next.front(), 0.0);
	EXPECT_NEAR(next[1], 0.5, 1e-12);
	EXPECT_DOUBLE_EQ(next.back(), 0.0);
}

TEST(FiniteMethods, GaussSeidelUsesUpdatedNeighbors) {
	FiniteMethodsFixture fx;
	fx.finite.setDiffusionModel(FiniteMethods::DiffusionModel::CONSTANT);
	fx.finite.setSolver(FiniteMethods::Solver::GAUSS_SEIDEL);
	fx.finite.setDiscretization(FiniteMethods::Discretization::FORWARD_EULER);
	fx.finite.setDiffusionCoefficient(1.0);
	fx.finite.setSpaceStep(1.0);
	fx.finite.setTimeStep(0.5); // alpha = 0.5
	std::vector<double> current{0.0, 1.0, 0.0, 0.0};
	std::vector<double> next;
	fx.finite._applyDiffusion(current, next);
	ASSERT_EQ(next.size(), current.size());
	EXPECT_DOUBLE_EQ(next[0], 0.0);
	EXPECT_NEAR(next[1], 0.0, 1e-12);
	EXPECT_NEAR(next[2], 0.0, 1e-12);
	EXPECT_DOUBLE_EQ(next[3], 0.0);
}
