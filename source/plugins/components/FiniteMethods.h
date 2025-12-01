/*
 * File:   FiniteMethods.h
 * Author: codex
 */

#ifndef FINITEMETHODS_H
#define FINITEMETHODS_H

#include <string>
#include <vector>
#include "../../kernel/simulator/ModelComponent.h"
#include "../data/Variable.h"

class FiniteMethods : public ModelComponent {
public:
	enum class DiffusionModel : int { NONE = 0, CONSTANT = 1, num_elements = 2 };
	enum class Discretization : int { FORWARD_EULER = 0, CRANK_NICOLSON = 1, num_elements = 2 };
	enum class Solver : int { EXPLICIT = 0, GAUSS_SEIDEL = 1, num_elements = 2 };

public:
	static std::string convertEnumToStr(DiffusionModel model);
	static std::string convertEnumToStr(Discretization discretization);
	static std::string convertEnumToStr(Solver solver);

public:
	FiniteMethods(Model* model, std::string name = "");
	virtual ~FiniteMethods() = default;

public:
	void setFieldVariable(Variable* field);
	Variable* getFieldVariable() const;
	void setTimeVariable(Variable* timeVariable);
	Variable* getTimeVariable() const;
	void setDiffusionModel(DiffusionModel model);
	DiffusionModel getDiffusionModel() const;
	void setDiscretization(Discretization discretization);
	Discretization getDiscretization() const;
	void setSolver(Solver solver);
	Solver getSolver() const;
	void setDiffusionCoefficient(double diffusionCoefficient);
	double getDiffusionCoefficient() const;
	void setSpaceStep(double spaceStep);
	double getSpaceStep() const;
	void setTimeStep(double timeStep);
	double getTimeStep() const;
	void setFilename(std::string filename);
	std::string getFilename() const;

public:
	virtual std::string show();

public:
	static PluginInformation* GetPluginInformation();
	static ModelComponent* LoadInstance(Model* model, PersistenceRecord *fields);
	static ModelDataDefinition* NewInstance(Model* model, std::string name = "");

protected:
	virtual void _onDispatchEvent(Entity* entity, unsigned int inputPortNumber);
	virtual bool _loadInstance(PersistenceRecord *fields);
	virtual void _saveInstance(PersistenceRecord *fields, bool saveDefaultValues);
	virtual bool _check(std::string& errorMessage);

private:
	bool _doStep();
	void _applyDiffusion(const std::vector<double>& current, std::vector<double>& next);

private:
	const struct DEFAULT_VALUES {
		const DiffusionModel diffusionModel = DiffusionModel::CONSTANT;
		const Discretization discretization = Discretization::FORWARD_EULER;
		const Solver solver = Solver::EXPLICIT;
		const double diffusionCoefficient = 0.5;
		const double spaceStep = 1.0;
		const double timeStep = 0.1;
		const std::string filename = "";
	} DEFAULT;

	Variable* _fieldVariable = nullptr;
	Variable* _timeVariable = nullptr;
	DiffusionModel _diffusionModel = DEFAULT.diffusionModel;
	Discretization _discretization = DEFAULT.discretization;
	Solver _solver = DEFAULT.solver;
	double _diffusionCoefficient = DEFAULT.diffusionCoefficient;
	double _spaceStep = DEFAULT.spaceStep;
	double _timeStep = DEFAULT.timeStep;
	std::string _filename = DEFAULT.filename;
};

#endif /* FINITEMETHODS_H */
